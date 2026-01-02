#include "graph.hpp"
#include "mfem.hpp"


using namespace mfem;

Graph::Graph(const FiniteElementSpace &fes) 
{
    int ne = fes.GetNE();
    graph = SparseMatrix(ne, ne);
    const Table &element_to_dof = fes.GetElementToDofTable();

    for(int i = 0; i < ne; i++)
    {
        node_to_element.Append(i);
        
        Array<int> idofs;
        element_to_dof.GetRow(i,idofs);
        
        for(int j = 0; j < ne; j++)
        {
            // An element should not be connected to itself.
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
};

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