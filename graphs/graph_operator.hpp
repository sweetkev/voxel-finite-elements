#pragma once

#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include "graph.hpp"
#include <unordered_map>

using namespace mfem;

class GraphOperator
{
public:
    GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_);

private:
    Array<int> A_ref;
    FiniteElementSpace reference_fes;

    // Matrix whose (i,j) entry denotes the local dof number on reference
    // element j which is identified with the local dof i on 
    // the central reference element.
    DenseMatrix local_to_neighbor_dof_map;
    
    SparseMatrix Q;
    Graph graph;

    /** Creates the table whose ith row represents a local dof number
     *  and jth column represents a reference cell in the reference mesh
     */
    void CreateDofMap();

    /** Returns the local dof number on neighbor element that corresponds 
     *  to the local dof number e_dof on element e. If the dof is not shared
     *  between elements, returns -1
    */
    int GetNeighborDof(int e, int e_dof, int neighbor);
};

void CreateReferenceMesh(Mesh &reference_mesh,int dim);