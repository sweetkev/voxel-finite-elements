#include "graph.hpp"
#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>

using namespace mfem;
using namespace std;

// Helper functions without reference to class members
namespace
{
/** Removes diagonal entries (representing self-connections) from graph */
void RemoveSelfConnections(Table &graph)
{
   int size = graph.Size();
   int *I = graph.GetI();
   int *J = graph.GetJ();

   // Count non-self entries (necessary for the case of an isolated node)
   int non_self_count = 0;
   for (int i = 0; i < size; ++i)
   {
      for (int j = I[i]; j < I[i+1]; ++j)
      {
         if (J[j] != i) { non_self_count++; }
      }
   }

   int *newI = new int[size + 1];
   int *newJ = new int[non_self_count];

   int idx = 0;
   for (int i = 0; i < size; ++i)
   {
      newI[i] = idx;
      for (int j = I[i]; j < I[i+1]; ++j)
      {
         if (J[j] == i) { continue; }
         newJ[idx] = J[j];
         idx++;
      }
   }
   newI[size] = idx;

   graph.SetIJ(newI, newJ, size);
}
}

// Class member functions

Graph::Graph(Table &graph_,
             unordered_map<Coord, vector<int>> &coord_to_node_,
             vector<Coord> &node_to_coord_)
   : graph(graph_),
     coord_to_node(coord_to_node_),
     node_to_coord(node_to_coord_) { }

Graph::Graph(const FiniteElementSpace &fes, const PixelImage &image_)
{
   // Element-to-element table is (element-to-dof) * (dof-to-element)
   const Table &element_to_dof = fes.GetElementToDofTable();
   unique_ptr<Table> dof_to_element(Transpose(element_to_dof));
   Mult(element_to_dof, *dof_to_element, graph);

   // Create maps between coordinates and vertices
   const int width = image_.Width();
   const int height = image_.Height();
   int e = 0;
   for (int j = 0; j < height; ++j)
   {
      for (int i = 0; i < width; ++i)
      {
         if (image_(i, j) != 0)
         {
            Coord coord(i, j);
            coord_to_node[coord].push_back(e);
            node_to_coord.push_back(coord);
            grid_cells.push_back(coord);
            ++e;
         }
      }
   }
}

Graph Graph::CoarsenGraph()
{
   // Create graph labeling and create node-cell maps over coarse graph
   vector<Coord> coarse_node_to_coord;
   unordered_map<Coord, vector<int>> coarse_coord_to_node;

   Array<int> graph_labeling = LabelGraph(coarse_coord_to_node,
                                          coarse_node_to_coord);

   // Create coarse graph from fine graph and its labeling
   Table coarse_graph = BuildCoarseGraph(graph_labeling,
                                         coarse_node_to_coord.size());

   return Graph(coarse_graph, coarse_coord_to_node, coarse_node_to_coord);
}

Array<int> Graph::LabelGraph(
   unordered_map<Coord, vector<int>> &coarse_coord_to_node,
   vector<Coord> &coarse_node_to_coord)
{
   // In each 2x2 block of cells that becomes a coarse cell, there are
   // coarse nodes for each disjoint connected set of fine nodes. For a
   // standard coarsening, this gives 1 coarse node for each 2x2 block of fine
   // cells.
   Array<int> graph_labeling(graph.Size());
   graph_labeling = -1;

   // Overlay coarse grid on top of fine grid
   unordered_map<Coord, vector<int>> coarse_coord_to_fine_node;
   vector<Coord> fine_node_to_coarse_coord(graph.Size());
   vector<Coord> coarse_grid_cells;

   for (Coord coord : grid_cells)
   {
      Coord coarse_coord(coord[0]/2, coord[1]/2);
      coarse_grid_cells.push_back(coarse_coord);
      for (int fine_node : coord_to_node[coord])
      {
         fine_node_to_coarse_coord[fine_node] = coarse_coord;
         coarse_coord_to_fine_node[coarse_coord].push_back(fine_node);
      }
   }

   int coarse_node_counter = 0;

   for (Coord coarse_coord : coarse_grid_cells)
   {
      const auto &fine_nodes = coarse_coord_to_fine_node[coarse_coord];
      const int n_fine_nodes = fine_nodes.size();
      vector<set<int>> connected_components(n_fine_nodes);
      unordered_map<int,int> fine_node_component_index(n_fine_nodes);

      auto merge = [&](int i, int j)
      {
         const int ci = fine_node_component_index[i];
         const int cj = fine_node_component_index[j];
         fine_node_component_index[j] = ci;
         connected_components[ci].merge(connected_components[cj]);
      };

      for (int i = 0; i < n_fine_nodes; ++i)
      {
         connected_components[i].emplace(fine_nodes[i]);
         fine_node_component_index[fine_nodes[i]] = i;
      }

      Array<int> row;
      for (int fine_node : fine_nodes)
      {
         graph.GetRow(fine_node, row);
         for (int neighbor_fine_node : row)
         {
            if (fine_node_to_coarse_coord[neighbor_fine_node] != coarse_coord)
            {
               continue;
            }
            merge(fine_node, neighbor_fine_node);
         }
      }

      vector<int> nonempty_component_index(n_fine_nodes);
      int nonempty_component_counter = 0;
      for (int i = 0; i < n_fine_nodes; ++i)
      {
         if (!connected_components[i].empty())
         {
            nonempty_component_index[i] = nonempty_component_counter;
            coarse_coord_to_node[coarse_coord].push_back(coarse_node_counter +
                                                         nonempty_component_counter);
            ++nonempty_component_counter;
         }
         else
         {
            nonempty_component_index[i] = -1;
         }
      }

      for (int fine_node : fine_nodes)
      {
         graph_labeling[fine_node] = coarse_node_counter +
                                     nonempty_component_index[fine_node_component_index[fine_node]];
      }

      coarse_node_counter += nonempty_component_counter;
   }

   coarse_node_to_coord.resize(coarse_node_counter);
   for (Coord coarse_coord : coarse_grid_cells)
   {
      for (int coarse_node : coarse_coord_to_node[coarse_coord])
      {
         coarse_node_to_coord[coarse_node] = coarse_coord;
      }
   }

   return graph_labeling;
}

Table Graph::BuildCoarseGraph(Array<int> &graph_labeling,
                              int coarse_ne)
{
   // Build connectivity matrix
   Table P(graph.Size(), 1);
   for (int i = 0; i < graph.Size(); ++i)
   {
      P.GetRow(i)[0] = graph_labeling[i];
   }

   unique_ptr<Table> Pt(Transpose(P));

   // Build coarse graph
   // G_coarse = P^T * G_fine * P
   Table coarse_graph;
   Table temp;

   Mult(*Pt, graph, temp);
   Mult(temp, P, coarse_graph);

   Table coarse_graph_fixed;
   coarse_graph_fixed.MakeI(coarse_graph.Size());
   Array<int> row;
   for (int i = 0; i < coarse_graph.Size(); ++i)
   {
      coarse_graph.GetRow(i, row);
      for (int j = 0; j < row.Size(); ++j)
      {
         if (row[j] != i) { coarse_graph_fixed.AddAColumnInRow(i); }
      }
   }
   coarse_graph_fixed.MakeJ();
   for (int i = 0; i < coarse_graph.Size(); ++i)
   {
      coarse_graph.GetRow(i, row);
      int k = 0;
      for (int j = 0; j < row.Size(); ++j)
      {
         if (row[j] != i)
         {
            coarse_graph_fixed.GetRow(i)[k] = row[j];
            ++k;
         }
      }
   }

   return coarse_graph_fixed;
}
