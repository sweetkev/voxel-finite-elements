#include "graph.hpp"
#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>


using namespace mfem;

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
        for(int i = 0; i < size; ++i)
        {
            for(int j = I[i]; j < I[i+1]; ++j)
            {
                if(J[j] != i) { non_self_count++; }
            }
        }

        int *newI = new int[size + 1];
        int *newJ = new int[non_self_count];
    
        int idx = 0;
        for(int i = 0; i < size; ++i)
        {
            newI[i] = idx;
            for(int j = I[i]; j < I[i+1]; ++j)
            {
                if(J[j] == i) { continue; }
                newJ[idx] = J[j];
                idx++;
            }
        }
        newI[size] = idx;

        graph.SetIJ(newI, newJ, size);
    }

    /** Returns true if the two arrays of DoFs share any DoFs */
    bool SharesDof(const Array<int> &idofs, const Array<int> &jdofs)
    {
        for(int i = 0; i < idofs.Size(); i++)
        {
            for(int j = 0; j < jdofs.Size(); j++)
            {
                if(idofs[i] == jdofs[j]) {
                    return true;
                }
            }
        }

        return false;
    }
}

// Class member functions

Graph::Graph(const FiniteElementSpace &fes, const PixelImage &image_) 
    : image(image_)
{
    int ne = fes.GetNE();
    graph = Table(ne, 8); // A fully surrounded element has at most 8 neighbors
    const Table &element_to_dof = fes.GetElementToDofTable();

    // Create graph
    for(int i = 0; i < ne; i++)
    {
        node_to_cell.Append(i);
        cell_to_node.push_back({i});

        Array<int> idofs;
        element_to_dof.GetRow(i,idofs);
        
        for(int j = 0; j < ne; j++)
        {
            // An cell should not be connected to itself
            if (i == j) { continue; }

            Array<int> jdofs;
            element_to_dof.GetRow(j,jdofs);
            
            if(SharesDof(idofs,jdofs)) {
                graph.Push(i,j);
                continue;
            }
        }
    }
    graph.Finalize();

    CreateCoordCellMaps(coord_to_cell, cell_to_coord, image);   
}

Graph Graph::CoarsenGraph() 
{
    // Generate coarsened geometry
    PixelImage coarse_image = image.Coarsen();

    // Create cell-coordinate maps over coarse graph
    std::unordered_map<Coord, int> coarse_coord_to_cell;
    std::vector<Coord> coarse_cell_to_coord;
    CreateCoordCellMaps(coarse_coord_to_cell, coarse_cell_to_coord,
        coarse_image);

    // Create graph labeling and create node-cell maps over coarse graph
    Array<int> coarse_node_to_cell;
    std::vector<std::vector<int>> 
        coarse_cell_to_node(coarse_cell_to_coord.size());
    Array<int> graph_labeling = LabelGraph(coarse_node_to_cell, 
        coarse_cell_to_node, coarse_coord_to_cell,
        coarse_cell_to_coord, coarse_image);

    // Create coarse graph from fine graph and its labeling
    Table coarse_graph = BuildCoarseGraph(graph_labeling,
        coarse_node_to_cell.Size());

    return Graph(coarse_graph,coarse_node_to_cell, coarse_cell_to_node,
        coarse_coord_to_cell, coarse_cell_to_coord, coarse_image);
}

void Graph::CreateCoordCellMaps(
    std::unordered_map<Coord, int> &coord_to_cell_,
    std::vector<Coord> &cell_to_coord_, const PixelImage &image_)
{
    // Create maps between coordinates and vertices
    int width = image_.Width(), height = image_.Height();
    cell_to_coord_.resize(width*height);
    int e = 0;
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
        if (image_(i, j) != 0)
        {
            Coord coord(i, j);
            coord_to_cell_[coord] = e;
            cell_to_coord_[e] = coord;
            e++;
        }
      }
   }
   cell_to_coord_.resize(e);   
}

/* Each child index corresponds to a fine cell. The sub_cell_position indicates
   the position of the fine cell in realtion to its parental coarse cell.
   sub_cell_position is lexicographic left-to-right, bottom-to-top, i.e.
   2 3
   0 1
*/
struct ChildIndex
{
    int fine_cell;
    int sub_cell_position;
};

Array<int> Graph::LabelGraph(
    Array<int> &coarse_node_to_cell, 
    std::vector<std::vector<int>> &coarse_cell_to_node,
    const std::unordered_map<Coord, int> &coarse_coord_to_cell,
    const std::vector<Coord> &coarse_cell_to_coord, 
    const PixelImage &coarse_image)
{
    // In each 2x2 block of cells that becomes a coarse cell, there are
    // coarse nodes for each disjoint connected set of fine nodes. For a
    // standard coarsening, this gives 1 coarse node for each 2x2 block of fine
    // cells.

    Array<int> graph_labeling;
    graph_labeling.SetSize(graph.Size(), -1);

    // int new_width = coarse_image.Width(), new_height = coarse_image.Height();
    int fine_ne = cell_to_coord.size();
    int coarse_ne = coarse_cell_to_coord.size();
    

    /* children and children_offsets utilize a CSR-like structure. For a
       coarse cell i, its children can be accessed by the indicies from
       children_offsets[i] up to (not including) children_offsets[i+1].
    */
    std::vector<ChildIndex> children(fine_ne);
    std::vector<int> children_offsets(coarse_ne + 1);
    int offset = 0;
    for (int i = 0; i < coarse_ne; ++i)
    {
        children_offsets[i] = offset;

        const Coord coarse_coord = coarse_cell_to_coord[i];
        Coord fine_coord(2*coarse_coord[0], 2*coarse_coord[1]);

        for (int jj = 0; jj < 2; ++jj)
        {
            fine_coord[1] += jj;
            for (int ii = 0; ii < 2; ++ii)
            {
                fine_coord[0] += ii;
                const int fine_idx = GetCellIndex(fine_coord);
                if (fine_idx >= 0)
                {
                    children[offset] = {fine_idx, ii + 2*jj};
                    ++offset;
                }
                fine_coord[0] -= ii;
            }
        fine_coord[1] -= jj;
        }
    }
    children_offsets.back() = offset;

    // Create labeling
    int label = 0;
    for(int i = 0; i < coarse_ne; ++i)
    {
        // Collect all fine nodes in the coarse cell's children
        std::vector<int> coarse_cell_nodes;
        for(int j = children_offsets[i]; j < children_offsets[i+1]; ++j)
        {
            int fine_cell = children[j].fine_cell;
            coarse_cell_nodes.insert(
                coarse_cell_nodes.end(),
                cell_to_node[fine_cell].begin(),
                cell_to_node[fine_cell].end());
        }

        // Add distinct label to each connected component among the fine nodes
        std::vector<bool> visited(coarse_cell_nodes.size(), false);
        for(int j = 0; j < coarse_cell_nodes.size(); ++j)
        {
            if(visited[j]) { continue; }

            graph_labeling[coarse_cell_nodes[j]] = label;
            visited[j] = true;
            for(int k = j+1; k < coarse_cell_nodes.size(); ++k)
            {
                if(visited[k]) { continue; }

                if(graph(coarse_cell_nodes[j],coarse_cell_nodes[k]) >= 0)
                {
                    graph_labeling[coarse_cell_nodes[k]] = label;
                    visited[k] = true;
                }
            }
            coarse_cell_to_node[i].push_back(label);
            coarse_node_to_cell.Append(i);
            label++;
        }
    }
    
    return graph_labeling;
}

Table Graph::BuildCoarseGraph(Array<int> &graph_labeling,
        int coarse_ne) 
{
    // Build connectivity matrix
    Table P(graph.Size(), coarse_ne);
    for(int i = 0; i < graph.Size(); ++i)
    {
        P.Push(i, graph_labeling[i]);
    }
    P.Finalize();

    // Build coarse graph
    // G_coarse = P^T * G_fine * P
    Table coarse_graph(coarse_ne, coarse_ne);
    Table Pt;
    Table temp;

    Transpose(P, Pt, coarse_ne);
    Pt.Finalize();
    
    Mult(Pt, graph, temp);
    Mult(temp, P, coarse_graph);

    // Remove self-connections in graph
    RemoveSelfConnections(coarse_graph);

    coarse_graph.Finalize();
    return coarse_graph;
}
