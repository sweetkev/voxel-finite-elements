#include "graph.hpp"
#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>


using namespace mfem;

Graph::Graph(const FiniteElementSpace &fes, const PixelImage &image_) 
    : image(image_)
{
    int ne = fes.GetNE();
    graph = SparseMatrix(ne, ne);
    const Table &element_to_dof = fes.GetElementToDofTable();

    // Create graph
    for(int i = 0; i < ne; i++)
    {
        node_to_element.Append(i);
        element_to_node.push_back({i});

        Array<int> idofs;
        element_to_dof.GetRow(i,idofs);
        
        for(int j = 0; j < ne; j++)
        {
            // An element should not be connected to itself
            if (i == j) { continue; }

            Array<int> jdofs;
            element_to_dof.GetRow(j,jdofs);
            
            if(SharesDof(idofs,jdofs)) {
                graph.Add(i,j,1.0);
                continue;
            }
        }
    }
    graph.Finalize();

    CreateCoordElementMaps(coord_to_element, element_to_coord, image);   
}

Graph Graph::CoarsenGraph() 
{
    // Generate coarsened geometry
    PixelImage coarse_image = image.Coarsen();

    // Create element-coordinate maps over coarse graph
    std::unordered_map<Coord, int> coarse_coord_to_element;
    std::vector<Coord> coarse_element_to_coord;
    CreateCoordElementMaps(coarse_coord_to_element, coarse_element_to_coord,
        coarse_image);

    // Create graph labeling and create node-element maps over coarse graph
    Array<int> coarse_node_to_element;
    std::vector<std::vector<int>> 
        coarse_element_to_node(coarse_element_to_coord.size());
    Array<int> graph_labeling = LabelGraph(coarse_node_to_element, 
        coarse_element_to_node, coarse_coord_to_element,
        coarse_element_to_coord, coarse_image);

    // TODO: Create coarse graph from fine graph and its labeling
    SparseMatrix coarse_graph;

    return Graph(coarse_graph,coarse_node_to_element, coarse_element_to_node,
        coarse_coord_to_element, coarse_element_to_coord, coarse_image);
}

bool Graph::SharesDof(const Array<int> &idofs, const Array<int> &jdofs)
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

void Graph::CreateCoordElementMaps(
    std::unordered_map<Coord, int> &coord_to_element_,
    std::vector<Coord> &element_to_coord_, const PixelImage &image_)
{
    // Create maps between coordinates and vertices
    int width = image_.Width(), height = image_.Height();
    element_to_coord_.resize(width*height);
    int e = 0;
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
        if (image_(i, j) != 0)
        {
            Coord coord(i, j);
            coord_to_element_[coord] = e;
            element_to_coord_[e] = coord;
            e++;
        }
      }
   }
   element_to_coord_.resize(e);   
}

/* Each child index corresponds to a fine element. The sub_element_position indicates
   the position of the fine element in realtion to its parental coarse element.
   sub_element_position is lexicographic left-to-right, bottom-to-top, i.e.
   2 3
   0 1
*/
struct ChildIndex
{
    int fine_element;
    int sub_element_position;
};

Array<int> Graph::LabelGraph(
    Array<int> &coarse_node_to_element, 
    std::vector<std::vector<int>> &coarse_element_to_node,
    const std::unordered_map<Coord, int> &coarse_coord_to_element,
    const std::vector<Coord> &coarse_element_to_coord, 
    const PixelImage &coarse_image)
{
    // In each 2x2 block of elements that becomes a coarse element, there are
    // coarse nodes for each disjoint connected set of fine nodes. For a
    // standard coarsening, this gives 1 coarse node for each 2x2 block of fine
    // elements.

    Array<int> graph_labeling;
    graph_labeling.SetSize(graph.Height(), -1);

    int new_width = coarse_image.Width(), new_height = coarse_image.Height();
    int fine_ne = element_to_coord.size();
    int coarse_ne = coarse_element_to_coord.size();
    

    /* children and children_offsets utilize a CSR-like structure. For a
       coarse element i, its children can be accessed by the indicies from
       children_offsets[i] up to (not including) children_offsets[i+1].
    */
    std::vector<ChildIndex> children(fine_ne);
    std::vector<int> children_offsets(coarse_ne + 1);
    int offset = 0;
    for (int i = 0; i < coarse_ne; ++i)
    {
        children_offsets[i] = offset;

        const Coord coarse_coord = coarse_element_to_coord[i];
        Coord fine_coord(2*coarse_coord[0], 2*coarse_coord[1]);

        for (int jj = 0; jj < 2; ++jj)
        {
            fine_coord[1] += jj;
            for (int ii = 0; ii < 2; ++ii)
            {
                fine_coord[0] += ii;
                const int fine_idx = GetElementIndex(fine_coord);
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
        // Collect all fine nodes in the coarse element's children
        std::vector<int> coarse_element_nodes;
        for(int j = children_offsets[i]; j < children_offsets[i+1]; ++j)
        {
            int fine_element = children[j].fine_element;
            coarse_element_nodes.insert(
                coarse_element_nodes.end(),
                element_to_node[fine_element].begin(),
                element_to_node[fine_element].end());
        }

        // Add distinct label to each connected component among the fine nodes
        std::vector<bool> visited(coarse_element_nodes.size(), false);
        for(int j = 0; j < coarse_element_nodes.size(); ++j)
        {
            if(visited[j]) { continue; }

            graph_labeling[coarse_element_nodes[j]] = label;
            visited[j] = true;
            for(int k = j+1; k < coarse_element_nodes.size(); ++k)
            {
                if(visited[k]) { continue; }

                if(graph(coarse_element_nodes[j],coarse_element_nodes[k]) > 0)
                {
                    graph_labeling[coarse_element_nodes[k]] = label;
                    visited[k] = true;
                }
            }
            coarse_element_to_node[i].push_back(label);
            coarse_node_to_element.Append(i);
            label++;
        }
    }
    
    return graph_labeling;
}


