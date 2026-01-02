#include "graph.hpp"
#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>


using namespace mfem;

Graph::Graph(const FiniteElementSpace &fes, const PixelImage &image) 
{
    int ne = fes.GetNE();
    graph = SparseMatrix(ne, ne);
    const Table &element_to_dof = fes.GetElementToDofTable();

    element_to_node.SetSize(ne,1);

    // Create graph
    for(int i = 0; i < ne; i++)
    {
        node_to_element.Append(i);
        element_to_node(i,0) = i;

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

    // Create maps between coordinates and vertices
    int width = image.Width(), height = image.Height();
    element_to_coord.resize(width*height);
    int e = 0;
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
        if (image(i, j) != 0)
        {
            Coord coord(i, j);
            coord_to_element[coord] = e;
            element_to_coord[e] = coord;
            e++;
        }
      }
   }    
}

Array<int> Graph::GetElementNodes(int i) {
    Array<int> element_nodes;
    element_to_node.GetRow(i,element_nodes);
    return element_nodes;
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