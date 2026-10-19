#include "voxel_graph.hpp"
#include "mfem.hpp"
#include "voxel_image.hpp"
#include <unordered_map>

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
}

VoxelGraph::VoxelGraph(const FiniteElementSpace &fes, const VoxelImage &image, const FiniteElementSpace &reference_fes_)
    : reference_fes(reference_fes_)
{
   // Element-to-element table is (element-to-dof) * (dof-to-element)
   const Table &element_to_dof = fes.GetElementToDofTable();
   unique_ptr<Table> dof_to_element(Transpose(element_to_dof));
   Table new_graph;
   Mult(element_to_dof, *dof_to_element, new_graph);
   graph = RemoveSelfConnections(new_graph);

   // Create maps between grid cells and elements
   element_to_cell = image.GetElementToCell();
   for (int e = 0; e < (int)element_to_cell.size(); ++e)
   {
      const CellIndex &cell = element_to_cell[e];
      vector<int> &elements = cell_to_elements[cell];
      if (elements.empty()) { grid_cells.push_back(cell); }
      elements.push_back(e);
   }

   CreateReferenceDofMap();
   CreateDofGroups();

   ReindexDofs(fes);
}

VoxelGraph VoxelGraph::CoarsenGraph()
{
    // Initialize coarse graph maps
    vector<CellIndex> coarse_grid_cells;
    unordered_map<CellIndex, vector<int>> coarse_cell_to_elements;
    vector<CellIndex> coarse_element_to_cell;

    // Initialize maps between coarse grid cells and fine elements
    vector<CellIndex> fine_element_to_coarse_cell(graph.Size());
    unordered_map<CellIndex, vector<int>> coarse_cell_to_fine_elements;

    // Find coarse grid cells
    for (CellIndex fine_cell : grid_cells)
    {
        CellIndex coarse_cell(fine_cell[0]/2, fine_cell[1]/2, fine_cell[2]/2);

        // Only add each coarse grid cell to list once
        const auto coarse_cell_visited = find(coarse_grid_cells.begin(),
                                              coarse_grid_cells.end(),
                                              coarse_cell);
        if (coarse_cell_visited == coarse_grid_cells.end())
        {
            coarse_grid_cells.push_back(coarse_cell);
        }

        // Create maps between the coarse grid cells and fine elements
        for (int fine_element : cell_to_elements[fine_cell])
        {
            fine_element_to_coarse_cell[fine_element] = coarse_cell;
            coarse_cell_to_fine_elements[coarse_cell].push_back(fine_element);
        }
    }

    // Label each fine element with its associated coarse element and create the
    // maps between coarse grid cells and coarse elements. This is done by
    // finding the connected components of the fine graph over each coarse grid
    // cell. Each connected component recieves a distinct label.
    int coarse_element_counter = 0;
    for (CellIndex coarse_cell : coarse_grid_cells)
    {
        const vector<int> &fine_elements = coarse_cell_to_fine_elements[coarse_cell];
        const int fine_ne = fine_elements.size();
        vector<set<int>> connected_components(fine_ne);
        unordered_map<int,int> fine_element_component_index(fine_ne);

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

        // Merge connected elements on coarse grid cell
        Array<int> row;
        for (int fine_element : fine_elements)
        {
            graph.GetRow(fine_element, row);
            for (int neighbor_fine_element : row)
            {
                if (fine_element_to_coarse_cell[neighbor_fine_element] != coarse_cell)
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
                coarse_cell_to_elements[coarse_cell].push_back(coarse_element_counter +
                                                         nonempty_component_counter);
                ++nonempty_component_counter;
            }
            else
            {
                nonempty_component_index[i] = -1;
            }
        }

        // Label graph nodes by their connected component
        fine_to_coarse_element_map.resize(graph.Size());
        for (int fine_element : fine_elements)
        {
            fine_to_coarse_element_map[fine_element] = coarse_element_counter +
                                         nonempty_component_index[fine_element_component_index[fine_element]];
        }

        coarse_element_counter += nonempty_component_counter;
    }

    // Create map between coarse elements and their grid cells
    coarse_element_to_cell.resize(coarse_element_counter);
    for (CellIndex coarse_cell : coarse_grid_cells)
    {
        for (int coarse_element : coarse_cell_to_elements[coarse_cell])
        {
            coarse_element_to_cell[coarse_element] = coarse_cell;
        }
   }

    // Create map between coarse elements and fine elements
    coarse_to_fine_elements_map.resize(coarse_element_counter);
    for (int fine_element = 0; fine_element < graph.Size(); ++fine_element)
    {
        int coarse_element = fine_to_coarse_element_map[fine_element];
        coarse_to_fine_elements_map[coarse_element].insert(fine_element);
    }

   // Build coarse graph according to coarse element labeling
   // Build connectivity matrix P
   Table P(graph.Size(), 1);
   for (int fine_element = 0; fine_element < graph.Size(); ++fine_element)
   {
      P.GetRow(fine_element)[0] = fine_to_coarse_element_map[fine_element];
   }

   unique_ptr<Table> Pt(Transpose(P));

   // Build coarse graph as G_coarse = P^T * G_fine * P
   Table coarse_graph_temp1;
   Table coarse_graph_temp2;

   Mult(*Pt, graph, coarse_graph_temp1);
   Mult(coarse_graph_temp1, P, coarse_graph_temp2);

   // Remove self-adjacencies
   Table coarse_graph = RemoveSelfConnections(coarse_graph_temp2);

   return VoxelGraph(coarse_graph, coarse_grid_cells, coarse_cell_to_elements,
                     coarse_element_to_cell, reference_fes);
}

VoxelGraph::VoxelGraph(Table &graph_, vector<CellIndex> &grid_cells_,
                         unordered_map<CellIndex, vector<int>> &cell_to_elements_,
                         vector<CellIndex> &element_to_cell_,
                         const FiniteElementSpace &reference_fes_)
    : graph(graph_), grid_cells(grid_cells_),
      cell_to_elements(cell_to_elements_), element_to_cell(element_to_cell_),
      reference_fes(reference_fes_)
{
    CreateReferenceDofMap();
    CreateDofGroups();
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

void VoxelGraph::CreateDofGroups()
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
        for (int e_dof = 0; e_dof < dofs_per_elem; ++e_dof)
        {
            Array<int> connected_elements;
            graph.GetRow(e, connected_elements);
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

    for (int i = 0; i < nbdofs; ++i)
    {
        if (!connected_components[i].empty())
        {
            dof_groups.push_back(connected_components[i]);
        }
    }

    // Build map from broken dof number to true dof number
    broken_to_true_dof.assign(nbdofs, -1);
    const int ntdofs = dof_groups.size();
    for (int tdof = 0; tdof < ntdofs; ++tdof)
    {
        const std::set<int> &dof_group = dof_groups[tdof];
        for (int broken_dof : dof_group)
        {
            broken_to_true_dof[broken_dof] = tdof;
        }
    }
}

int VoxelGraph::GetNeighborDof(int e, int e_dof, int neighbor)
{
    CellIndex e_cell = element_to_cell[e];
    CellIndex neighbor_cell = element_to_cell[neighbor];

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
    int ntdofs = dof_groups.size();
    int ne = graph.Size();
    int nbdofs = dofs_per_elem * ne;
    vector<set<int>> new_dof_groups(ntdofs);
    vector<int> voxel_to_mfem(ntdofs);
    for (int mfem_dof = 0; mfem_dof < ntdofs; ++mfem_dof)
    {
        int e = fes.GetElementForDof(mfem_dof);
        int local_mfem_index = fes.GetLocalDofForDof(mfem_dof);
        int broken_dof = e * dofs_per_elem + local_mfem_index;
        int voxel_dof = broken_to_true_dof[broken_dof];
        voxel_to_mfem[voxel_dof] = mfem_dof;

        new_dof_groups[mfem_dof] = dof_groups[voxel_dof];
    }
    dof_groups = new_dof_groups;

    for (int bdof = 0; bdof < nbdofs; ++bdof)
    {
        int voxel_dof = broken_to_true_dof[bdof];
        int mfem_dof = voxel_to_mfem[voxel_dof];
        broken_to_true_dof[bdof] = mfem_dof;
    }
}
