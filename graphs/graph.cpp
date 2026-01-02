#include "graph.hpp"
#include "mfem.hpp"


using namespace mfem;

Graph::Graph(const FiniteElementSpace &fes) 
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