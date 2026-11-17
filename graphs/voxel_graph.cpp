#include "voxel_graph.hpp"
#include "mfem.hpp"
#include "voxel_image.hpp"
#include <unordered_map>
#include <algorithm>

using namespace mfem;
using namespace std;

namespace
{
/** Removes diagonal entries (representing self-connections) from graph */
Table RemoveSelfConnections(Table &graph)
{
   Table graph_fixed;
   graph_fixed.MakeI(graph.Size());
   Array<int> row;
   for (int i = 0; i < graph.Size(); ++i)
   {
      graph.GetRow(i, row);
      for (int j = 0; j < row.Size(); ++j)
      {
         if (row[j] != i) { graph_fixed.AddAColumnInRow(i); }
      }
   }
   graph_fixed.MakeJ();
   for (int i = 0; i < graph.Size(); ++i)
   {
      graph.GetRow(i, row);
      int k = 0;
      for (int j = 0; j < row.Size(); ++j)
      {
         if (row[j] != i)
         {
            graph_fixed.GetRow(i)[k] = row[j];
            ++k;
         }
      }
   }
   return graph_fixed;
}

// Returns offsets with offsets[r] the global number of the first of the n
// elements on rank r, and offsets[nranks] the global count.
vector<int> GlobalOffsets(MPI_Comm comm, int n)
{
   int nranks;
   MPI_Comm_size(comm, &nranks);
   vector<int> counts(nranks);
   MPI_Allgather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
   vector<int> offsets(nranks + 1, 0);
   long long total = 0;
   for (int r = 0; r < nranks; ++r)
   {
      total += counts[r];
      offsets[r + 1] = total;
   }
   return offsets;
}

// Returns the rank holding the graph element with global number i.
int OffsetRank(const vector<int> &offsets, int i)
{
   return upper_bound(offsets.begin(), offsets.end(), i) - offsets.begin() - 1;
}

// Converts offsets from GlobalOffsets into a hypre partitioning array:
// [offsets[rank], offsets[rank + 1]] if HYPRE_AssumedPartitionCheck(), and all
// of offsets otherwise.
Array<HYPRE_BigInt> HyprePartitioning(const vector<int> &offsets, int rank)
{
   Array<HYPRE_BigInt> partitioning;
   if (HYPRE_AssumedPartitionCheck())
   {
      partitioning.Append(offsets[rank]);
      partitioning.Append(offsets[rank + 1]);
   }
   else
   {
      for (int offset : offsets) { partitioning.Append(offset); }
   }
   return partitioning;
}
}

VoxelGraph::VoxelGraph(const ParFiniteElementSpace &pfes,
                       vector<CellIndex> element_to_cell_,
                       const FiniteElementSpace &reference_fes_)
    : comm(pfes.GetComm()), rank(pfes.GetMyRank()),
      element_to_cell(std::move(element_to_cell_)),
      reference_fes(reference_fes_)
{
   // Element-to-element table is (element-to-dof) * (dof-to-element)
   const Table &element_to_dof = pfes.GetElementToDofTable();
   unique_ptr<Table> dof_to_element(Transpose(element_to_dof));
   Table new_graph;
   Mult(element_to_dof, *dof_to_element, new_graph);
   graph = RemoveSelfConnections(new_graph);

   element_offsets = GlobalOffsets(comm, graph.Size());
   FindOffRankNeighbors(pfes);

   // Create maps between grid cells and elements
   for (int e = 0; e < (int)element_to_cell.size(); ++e)
   {
      const CellIndex &cell = element_to_cell[e];
      vector<int> &elements = cell_to_elements[cell];
      if (elements.empty()) { grid_cells.push_back(cell); }
      elements.push_back(e);
   }

   CreateReferenceDofMap();
   CreateBrokenToLocalDofMap();

   ReindexDofs(pfes);
   dof_truedof = pfes.Dof_TrueDof_Matrix();
}

int VoxelGraph::ElementRank(int global_e) const
{
   return OffsetRank(element_offsets, global_e);
}

void VoxelGraph::FindOffRankNeighbors(const ParFiniteElementSpace &pfes)
{
   ParMesh &pmesh = *pfes.GetParMesh();
   const GroupTopology &gtopo = pmesh.gtopo;
   const int nranks = element_offsets.size() - 1;
   unique_ptr<Table> vertex_to_element(pmesh.GetVertexToElementTable());

   // Identify each shared vertex across ranks by the global true dof number
   // of its vertex dof.
   auto global_vertex = [&](int v)
   {
      Array<int> dofs;
      pfes.GetVertexDofs(v, dofs);
      const HYPRE_BigInt gv = pfes.GetGlobalTDofNumber(dofs[0]);
      return (int)gv;
   };

   // Send [global vertex, global element] for every local element at every
   // shared vertex to every other rank sharing that vertex. Group 0 is the
   // local group.
   vector<vector<int>> send(nranks);
   unordered_map<int, int> shared_vertex; // global vertex -> local vertex
   for (int group = 1; group < pmesh.GetNGroups(); ++group)
   {
      const int *group_ranks = gtopo.GetGroup(group);
      const int group_size = gtopo.GetGroupSize(group);
      for (int i = 0; i < pmesh.GroupNVertices(group); ++i)
      {
         const int v = pmesh.GroupVertex(group, i);
         const int gv = global_vertex(v);
         shared_vertex[gv] = v;
         for (int j = 0; j < group_size; ++j)
         {
            const int r = gtopo.GetNeighborRank(group_ranks[j]);
            if (r == rank) { continue; }
            const int *row = vertex_to_element->GetRow(v);
            for (int k = 0; k < vertex_to_element->RowSize(v); ++k)
            {
               send[r].insert(send[r].end(),
                              {gv, GlobalElementNumber(row[k])});
            }
         }
      }
   }

   vector<int> recv_counts;
   const vector<int> recv = Exchange(comm, send, recv_counts);

   // Every local element at a received vertex neighbors the remote element.
   Array<Connection> edges;
   for (size_t pos = 0; pos < recv.size(); pos += 2)
   {
      const int v = shared_vertex.at(recv[pos]);
      const int remote_e = recv[pos + 1];
      const int *row = vertex_to_element->GetRow(v);
      for (int k = 0; k < vertex_to_element->RowSize(v); ++k)
      {
         edges.Append(Connection(row[k], remote_e));
      }
   }
   edges.Sort();
   edges.Unique();
   off_rank_neighbors.MakeFromList(graph.Size(), edges);
}

VoxelGraph VoxelGraph::CoarsenGraph(const CellOwnership &ownership,
                                    int coarse_level)
{
   const int nranks = element_offsets.size() - 1;
   const int fine_level = coarse_level + 1;
   const int ne = graph.Size();

   // Send every fine element to the owner of its coarse grid cell as
   // [global element, i, j, k, n, neighbor_1, ..., neighbor_n], with its
   // fine grid cell (i, j, k) and its neighbors in global numbering.
   vector<vector<int>> to_owner(nranks);
   for (int e = 0; e < ne; ++e)
   {
      const CellIndex &fine_cell = element_to_cell[e];
      const CellIndex coarse_cell =
         AncestorCell(fine_cell, fine_level, coarse_level);
      vector<int> &buf = to_owner[ownership.GetOwner(coarse_level, coarse_cell)];
      const int *row = graph.GetRow(e);
      const int nlocal = graph.RowSize(e);
      const int *off_rank_row = off_rank_neighbors.GetRow(e);
      const int noff_rank = off_rank_neighbors.RowSize(e);
      buf.insert(buf.end(), {GlobalElementNumber(e), fine_cell[0],
                             fine_cell[1], fine_cell[2], nlocal + noff_rank});
      for (int k = 0; k < nlocal; ++k)
      {
         buf.push_back(GlobalElementNumber(row[k]));
      }
      buf.insert(buf.end(), off_rank_row, off_rank_row + noff_rank);
   }
   vector<int> recv_counts;
   const vector<int> recv = Exchange(comm, to_owner, recv_counts);
   to_owner.clear();

   // The fine elements received by this owner are numbered 0, ..., nfine - 1
   // in the order received. They form the owned fine graph, whose rows are
   // the owned fine elements and whose columns are the owned fine elements
   // followed by "ghost" neighbors: fine elements in cells of other owners,
   // numbered nfine, nfine + 1, ...
   vector<int> fine_global;           
   vector<CellIndex> fine_cells;
   unordered_map<int, int> fine_index; // global number -> owned fine element
   for (size_t pos = 0; pos < recv.size(); pos += 5 + recv[pos + 4])
   {
      fine_index[recv[pos]] = fine_global.size();
      fine_global.push_back(recv[pos]);
      fine_cells.emplace_back(recv[pos + 1], recv[pos + 2], recv[pos + 3]);
   }
   const int nfine = fine_global.size();

   vector<int> ghost_global; // global number of each ghost neighbor
   Array<Connection> fine_edges;
   for (size_t pos = 0, i = 0; pos < recv.size(); pos += 5 + recv[pos + 4], ++i)
   {
      for (int k = 0; k < recv[pos + 4]; ++k)
      {
         const int neighbor = recv[pos + 5 + k];
         auto it = fine_index.find(neighbor);
         if (it == fine_index.end())
         {
            it = fine_index.emplace(neighbor, nfine + ghost_global.size()).first;
            ghost_global.push_back(neighbor);
         }
         fine_edges.Append(Connection(i, it->second));
      }
   }
   fine_edges.Sort();
   const Table fine_graph(nfine, fine_edges);
   const int nghost = ghost_global.size();

   // Initialize coarse graph maps
   vector<CellIndex> coarse_grid_cells;
   unordered_map<CellIndex, vector<int>> coarse_cell_to_elements;
   vector<CellIndex> coarse_element_to_cell;

   // Initialize maps between coarse grid cells and owned fine elements
   vector<CellIndex> fine_element_to_coarse_cell(nfine);
   unordered_map<CellIndex, vector<int>> coarse_cell_to_fine_elements;

   // Find coarse grid cells
   for (int fine_element = 0; fine_element < nfine; ++fine_element)
   {
      const CellIndex coarse_cell =
         AncestorCell(fine_cells[fine_element], fine_level, coarse_level);

      // Only add each coarse grid cell to list once
      vector<int> &coarse_cell_fine_elements =
         coarse_cell_to_fine_elements[coarse_cell];
      if (coarse_cell_fine_elements.empty())
      {
         coarse_grid_cells.push_back(coarse_cell);
      }

      // Create maps between the coarse grid cells and fine elements
      fine_element_to_coarse_cell[fine_element] = coarse_cell;
      coarse_cell_fine_elements.push_back(fine_element);
   }

   // Label each owned fine element with its associated coarse element and
   // create the maps between coarse grid cells and coarse elements. This is
   // done by finding the connected components of the fine graph over each
   // coarse grid cell. Each connected component receives a distinct label.
   vector<int> fine_to_coarse(nfine); // local coarse numbering
   int coarse_element_counter = 0;
   for (const CellIndex &coarse_cell : coarse_grid_cells)
   {
      const vector<int> &fine_elements = coarse_cell_to_fine_elements[coarse_cell];
      const int fine_ne = fine_elements.size();
      vector<set<int>> connected_components(fine_ne);
      unordered_map<int, int> fine_element_component_index(fine_ne);

      auto merge = [&](int i, int j)
      {
         const int ci = fine_element_component_index[i];
         const int cj = fine_element_component_index[j];
         if (ci == cj) { return; }
         for (int member : connected_components[cj])
         {
            fine_element_component_index[member] = ci;
         }
         connected_components[ci].merge(connected_components[cj]);
         connected_components[cj].clear();
      };

      // Initialize each element with its own label
      for (int i = 0; i < fine_ne; ++i)
      {
         connected_components[i].emplace(fine_elements[i]);
         fine_element_component_index[fine_elements[i]] = i;
      }

      // Merge connected elements on coarse grid cell. Ghost neighbors lie in
      // other owners' cells, so they are never in this cell.
      Array<int> row;
      for (int fine_element : fine_elements)
      {
         fine_graph.GetRow(fine_element, row);
         for (int neighbor_fine_element : row)
         {
            if (neighbor_fine_element >= nfine ||
                fine_element_to_coarse_cell[neighbor_fine_element] != coarse_cell)
            {
               continue;
            }
            merge(fine_element, neighbor_fine_element);
         }
      }

      // Keep only the connected components which are nonempty
      vector<int> nonempty_component_index(fine_ne);
      int nonempty_component_counter = 0;
      for (int i = 0; i < fine_ne; ++i)
      {
         if (!connected_components[i].empty())
         {
            nonempty_component_index[i] = nonempty_component_counter;
            coarse_cell_to_elements[coarse_cell].push_back(
               coarse_element_counter + nonempty_component_counter);
            ++nonempty_component_counter;
         }
         else
         {
            nonempty_component_index[i] = -1;
         }
      }

      // Label graph nodes by their connected component
      for (int fine_element : fine_elements)
      {
         fine_to_coarse[fine_element] = coarse_element_counter +
            nonempty_component_index[fine_element_component_index[fine_element]];
      }

      coarse_element_counter += nonempty_component_counter;
   }
   const int coarse_ne = coarse_element_counter;

   // Create map between coarse elements and their grid cells
   coarse_element_to_cell.resize(coarse_ne);
   for (const CellIndex &coarse_cell : coarse_grid_cells)
   {
      for (int coarse_element : coarse_cell_to_elements[coarse_cell])
      {
         coarse_element_to_cell[coarse_element] = coarse_cell;
      }
   }

   // Create map between coarse elements and (global) fine elements
   coarse_to_fine_elements_map.assign(coarse_ne, {});
   for (int fine_element = 0; fine_element < nfine; ++fine_element)
   {
      coarse_to_fine_elements_map[fine_to_coarse[fine_element]].insert(
         fine_global[fine_element]);
   }

   // Global numbering of the coarse elements
   vector<int> coarse_offsets = GlobalOffsets(comm, coarse_ne);
   const int coarse_offset = coarse_offsets[rank];

   // Tell the rank holding each fine element its global coarse element as
   // [global fine element, global coarse element].
   vector<vector<int>> assignments(nranks);
   for (int fine_element = 0; fine_element < nfine; ++fine_element)
   {
      const int global_e = fine_global[fine_element];
      vector<int> &buf = assignments[ElementRank(global_e)];
      buf.insert(buf.end(),
                 {global_e, coarse_offset + fine_to_coarse[fine_element]});
   }
   const vector<int> recv_assignments = Exchange(comm, assignments, recv_counts);
   assignments.clear();

   fine_to_coarse_element_map.assign(ne, -1);
   for (size_t pos = 0; pos < recv_assignments.size(); pos += 2)
   {
      fine_to_coarse_element_map[recv_assignments[pos] - element_offsets[rank]] =
         recv_assignments[pos + 1];
   }

   // Ask the ranks holding the ghost neighbors for their global coarse
   // elements. Replies come back in the order of the queries.
   vector<vector<int>> queries(nranks);
   for (int neighbor : ghost_global)
   {
      queries[ElementRank(neighbor)].push_back(neighbor);
   }
   const vector<int> recv_queries = Exchange(comm, queries, recv_counts);

   vector<vector<int>> answers(nranks);
   for (int src = 0, pos = 0; src < nranks; ++src)
   {
      for (const int end = pos + recv_counts[src]; pos < end; ++pos)
      {
         const int e = recv_queries[pos] - element_offsets[rank];
         answers[src].push_back(fine_to_coarse_element_map[e]);
      }
   }
   const vector<int> recv_answers = Exchange(comm, answers, recv_counts);
   answers.clear();

   unordered_map<int, int> ghost_coarse; // global ghost -> global coarse
   for (int src = 0, pos = 0; src < nranks; ++src)
   {
      for (int neighbor : queries[src])
      {
         ghost_coarse[neighbor] = recv_answers[pos++];
      }
   }
   queries.clear();

   // Build coarse graph according to coarse element labeling. Columns of the
   // coarse graph are the owned coarse elements followed by the distinct
   // coarse elements of the ghost neighbors, which live on other ranks.
   vector<int> off_rank_coarse;                 // column - coarse_ne -> global
   unordered_map<int, int> off_rank_coarse_col; // global -> column
   Table P(nfine + nghost, 1);
   for (int fine_element = 0; fine_element < nfine; ++fine_element)
   {
      P.GetRow(fine_element)[0] = fine_to_coarse[fine_element];
   }
   for (int g = 0; g < nghost; ++g)
   {
      const int coarse_global = ghost_coarse.at(ghost_global[g]);
      auto it = off_rank_coarse_col.find(coarse_global);
      if (it == off_rank_coarse_col.end())
      {
         it = off_rank_coarse_col.emplace(coarse_global,
                                          coarse_ne + off_rank_coarse.size())
                                 .first;
         off_rank_coarse.push_back(coarse_global);
      }
      P.GetRow(nfine + g)[0] = it->second;
   }

   // Pt only needs the owned rows of P, which map to owned coarse elements.
   Table P_owned(nfine, 1);
   for (int fine_element = 0; fine_element < nfine; ++fine_element)
   {
      P_owned.GetRow(fine_element)[0] = fine_to_coarse[fine_element];
   }
   Table Pt;
   Transpose(P_owned, Pt, coarse_ne);

   // Build coarse graph as G_coarse = P^T * G_fine * P
   Table coarse_graph_temp1;
   Table coarse_graph_temp2;

   Mult(Pt, fine_graph, coarse_graph_temp1);
   Mult(coarse_graph_temp1, P, coarse_graph_temp2);

   // Remove self-adjacencies
   const Table coarse_graph_full = RemoveSelfConnections(coarse_graph_temp2);

   // Split the columns into local edges and off-rank neighbors (global
   // numbering).
   Array<Connection> local_edges, off_rank_edges;
   Array<int> row;
   for (int coarse_e = 0; coarse_e < coarse_ne; ++coarse_e)
   {
      coarse_graph_full.GetRow(coarse_e, row);
      for (int col : row)
      {
         if (col < coarse_ne)
         {
            local_edges.Append(Connection(coarse_e, col));
         }
         else
         {
            off_rank_edges.Append(Connection(coarse_e,
                                             off_rank_coarse[col - coarse_ne]));
         }
      }
   }
   local_edges.Sort();
   off_rank_edges.Sort();

   return VoxelGraph(comm, std::move(coarse_offsets),
                     Table(coarse_ne, local_edges),
                     Table(coarse_ne, off_rank_edges),
                     std::move(coarse_grid_cells),
                     std::move(coarse_cell_to_elements),
                     std::move(coarse_element_to_cell), reference_fes);
}

VoxelGraph::VoxelGraph(MPI_Comm comm_, vector<int> element_offsets_,
                       Table graph_, Table off_rank_neighbors_,
                       vector<CellIndex> grid_cells_,
                       unordered_map<CellIndex, vector<int>> cell_to_elements_,
                       vector<CellIndex> element_to_cell_,
                       const FiniteElementSpace &reference_fes_)
    : comm(comm_),
      element_offsets(std::move(element_offsets_)),
      grid_cells(std::move(grid_cells_)),
      cell_to_elements(std::move(cell_to_elements_)),
      element_to_cell(std::move(element_to_cell_)),
      reference_fes(reference_fes_)
{
    MPI_Comm_rank(comm, &rank);
    graph.Swap(graph_);
    off_rank_neighbors.Swap(off_rank_neighbors_);
    CreateReferenceDofMap();
    CreateBrokenToLocalDofMap();
    CreateDofTrueDofMatrix();
}

void VoxelGraph::CreateReferenceDofMap()
{
    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int d = reference_fes.GetMesh()->Dimension();

    local_to_neighbor_dof_map.SetSize(dofs_per_elem, pow(3,d));
    local_to_neighbor_dof_map = -1;

    Array<int> ref_dofs;
    int central_element = 0.5 * (pow(3,d)-1);
    reference_fes.GetElementDofs(central_element, ref_dofs);

    // For each local dof of reference element, find the neighboring elements
    // who share the dof, and list the corresponding local dof number on
    // neighbor.
    for (int i = 0; i < dofs_per_elem; ++i)
    {
        int dof = ref_dofs[i];
        for (int e = 0; e < pow(3,d); ++e)
        {
            if (e == central_element) { continue; }
            Array<int> neighbor_dofs;
            reference_fes.GetElementDofs(e, neighbor_dofs);
            for (int j = 0; j < neighbor_dofs.Size(); ++j)
            {
                if (neighbor_dofs[j] == dof)
                {
                    local_to_neighbor_dof_map(i,e) = j;
                    break;
                }
            }
        }
    }
}

void VoxelGraph::CreateBrokenToLocalDofMap()
{
    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int ne = graph.Size();
    const int nbdofs = dofs_per_elem * ne;
    Array<int> dof_labeling;

    std::vector<int> element_component_index(nbdofs);
    std::vector<std::set<int>> connected_components(nbdofs);

    // Initialize E_i = { i } and component index
    for (int i = 0; i < nbdofs; ++i)
    {
        connected_components[i].emplace(i);
        element_component_index[i] = i;
    }

    // Method to merge connected component j into connected component i
    auto merge = [&](int i, int j)
    {
        const int ci = element_component_index[i];
        const int cj = element_component_index[j];
        if (ci == cj) { return; }
        for (int member : connected_components[cj])
        {
            element_component_index[member] = ci;
        }
        connected_components[ci].merge(connected_components[cj]);
        connected_components[cj].clear();
    };

    for (int e = 0; e < ne; ++e)
    {
        Array<int> connected_elements;
        graph.GetRow(e, connected_elements);

        for (int e_dof = 0; e_dof < dofs_per_elem; ++e_dof)
        {
            for (int neighbor : connected_elements)
            {
                int neighbor_dof = GetNeighborDof(e, e_dof, neighbor);
                if (neighbor_dof != -1)
                {
                    int e_index = dofs_per_elem * e + e_dof;
                    int neighbor_index = dofs_per_elem * neighbor
                                         + neighbor_dof;
                    int ci = element_component_index[e_index];
                    int cj = element_component_index[neighbor_index];
                    if (ci != cj)
                    {
                        merge(e_index, neighbor_index);
                    }
                }
            }
        }
    }

    // Build map from broken dof number to local dof number. Each nonempty
    // connected component becomes one local dof.
    broken_to_local_dof.assign(nbdofs, -1);
    int dof = 0;
    for (int i = 0; i < nbdofs; ++i)
    {
        if (connected_components[i].empty()) { continue; }
        for (int broken_dof : connected_components[i])
        {
            broken_to_local_dof[broken_dof] = dof;
        }
        ++dof;
    }
    num_local_dofs = dof;
}

void VoxelGraph::CreateDofTrueDofMatrix()
{
    const int nranks = element_offsets.size() - 1;
    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int ndofs = num_local_dofs;

    // 1. Label each local dof with its global number.
    const vector<int> dof_offsets = GlobalOffsets(comm, ndofs);
    vector<int> label(ndofs);
    for (int dof = 0; dof < ndofs; ++dof)
    {
        label[dof] = dof_offsets[rank] + dof;
    }

    // 2. Lower each label to the smallest label of the matching dofs on
    //    off-rank neighbors until no label changes anywhere.
    // Sends [neighbor global element, i, j, k, labels of the
    //        element's dofs]
    const int message_size = 4 + dofs_per_elem;
    bool changed = true;
    while (changed)
    {
        vector<vector<int>> send(nranks);
        Array<int> neighbors;
        for (int e = 0; e < Size(); ++e)
        {
            const CellIndex &cell = element_to_cell[e];
            off_rank_neighbors.GetRow(e, neighbors);
            for (int neighbor : neighbors)
            {
                vector<int> &buf = send[ElementRank(neighbor)];
                buf.insert(buf.end(), {neighbor, cell[0], cell[1], cell[2]});
                for (int i = 0; i < dofs_per_elem; ++i)
                {
                    buf.push_back(label[broken_to_local_dof[e * dofs_per_elem + i]]);
                }
            }
        }
        vector<int> recv_counts;
        const vector<int> recv = Exchange(comm, send, recv_counts);

        int local_changed = 0;
        for (size_t pos = 0; pos < recv.size(); pos += message_size)
        {
            const int e = recv[pos] - element_offsets[rank];
            const CellIndex neighbor_cell(recv[pos + 1], recv[pos + 2],
                                          recv[pos + 3]);
            const int *neighbor_labels = &recv[pos + 4];
            for (int e_dof = 0; e_dof < dofs_per_elem; ++e_dof)
            {
                const int neighbor_dof =
                    GetNeighborDof(element_to_cell[e], e_dof, neighbor_cell);
                if (neighbor_dof == -1) { continue; }
                int &l = label[broken_to_local_dof[e * dofs_per_elem + e_dof]];
                if (neighbor_labels[neighbor_dof] < l)
                {
                    l = neighbor_labels[neighbor_dof];
                    local_changed = 1;
                }
            }
        }
        int global_changed = 0;
        MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_LOR, comm);
        changed = global_changed;
    }

    // 3. A local dof still labeled with its own global number holds the
    //    smallest label among the copies of its true dof, so this rank owns
    //    that true dof.
    auto owned = [&](int dof) { return label[dof] == dof_offsets[rank] + dof; };

    vector<int> dof_to_truedof(ndofs, -1);
    int ntruedofs = 0;
    for (int dof = 0; dof < ndofs; ++dof)
    {
        if (owned(dof)) { dof_to_truedof[dof] = ntruedofs++; }
    }
    const vector<int> truedof_offsets = GlobalOffsets(comm, ntruedofs);
    for (int dof = 0; dof < ndofs; ++dof)
    {
        if (owned(dof)) { dof_to_truedof[dof] += truedof_offsets[rank]; }
    }

    // Ask the rank holding the smallest label of each remaining dof for its
    // true dof number. The answers from each rank come back in the order of
    // the queries sent to it.
    // Sends [label]
    vector<vector<int>> queries(nranks);
    for (int dof = 0; dof < ndofs; ++dof)
    {
        if (!owned(dof))
        {
            queries[OffsetRank(dof_offsets, label[dof])].push_back(label[dof]);
        }
    }
    vector<int> query_counts;
    const vector<int> recv_queries = Exchange(comm, queries, query_counts);

    // Sends [true dof]
    vector<vector<int>> answers(nranks);
    {
        int pos = 0;
        for (int src = 0; src < nranks; ++src)
        {
            for (int k = 0; k < query_counts[src]; ++k, ++pos)
            {
                const int dof = recv_queries[pos] - dof_offsets[rank];
                answers[src].push_back(dof_to_truedof[dof]);
            }
        }
    }
    vector<int> answer_counts;
    const vector<int> recv_answers = Exchange(comm, answers, answer_counts);

    vector<int> next_answer(nranks, 0);
    for (int r = 1; r < nranks; ++r)
    {
        next_answer[r] = next_answer[r - 1] + answer_counts[r - 1];
    }
    for (int dof = 0; dof < ndofs; ++dof)
    {
        if (!owned(dof))
        {
            const int r = OffsetRank(dof_offsets, label[dof]);
            dof_to_truedof[dof] = recv_answers[next_answer[r]++];
        }
    }

    // Each row (local dof) has a single 1 in the column of its true dof.
    Array<int> I(ndofs + 1);
    Array<HYPRE_BigInt> J(ndofs);
    Vector data(ndofs);
    for (int dof = 0; dof < ndofs; ++dof)
    {
        I[dof] = dof;
        J[dof] = dof_to_truedof[dof];
    }
    I[ndofs] = ndofs;
    data = 1.0;

    const Array<HYPRE_BigInt> rows = HyprePartitioning(dof_offsets, rank);
    const Array<HYPRE_BigInt> cols = HyprePartitioning(truedof_offsets, rank);
    owned_dof_truedof = make_unique<HypreParMatrix>(
        comm, ndofs, dof_offsets[nranks], truedof_offsets[nranks],
        I.GetData(), J.GetData(), data.GetData(), rows.GetData(),
        cols.GetData());
    dof_truedof = owned_dof_truedof.get();
}

int VoxelGraph::GetNeighborDof(const CellIndex &e_cell, int e_dof,
                               const CellIndex &neighbor_cell) const
{
    int d = reference_fes.GetMesh()->Dimension();
    int neighbor_ref_element = 0;
    for (int i = 0; i < d; ++i)
    {
        neighbor_ref_element += (neighbor_cell[i] - e_cell[i] + 1) * pow(3, i);
    }
    return local_to_neighbor_dof_map(e_dof, neighbor_ref_element);
}


void VoxelGraph::ReindexDofs(const FiniteElementSpace &fes)
{
    int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    int nldofs = num_local_dofs;
    int ne = graph.Size();
    int nbdofs = dofs_per_elem * ne;
    vector<int> voxel_to_mfem(nldofs);
    for (int mfem_dof = 0; mfem_dof < nldofs; ++mfem_dof)
    {
        int e = fes.GetElementForDof(mfem_dof);
        int local_mfem_index = fes.GetLocalDofForDof(mfem_dof);
        int broken_dof = e * dofs_per_elem + local_mfem_index;
        int voxel_dof = broken_to_local_dof[broken_dof];
        voxel_to_mfem[voxel_dof] = mfem_dof;
    }

    for (int bdof = 0; bdof < nbdofs; ++bdof)
    {
        int voxel_dof = broken_to_local_dof[bdof];
        int mfem_dof = voxel_to_mfem[voxel_dof];
        broken_to_local_dof[bdof] = mfem_dof;
    }
}
