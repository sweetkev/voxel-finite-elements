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
    Table dof_to_reference_element_map;
    SparseMatrix Q;
    Graph graph;

    /** Creates the table whose ith row represents a local dof number
     *  and jth column represents a reference cell in the reference mesh
     */
    void CreateDofToElementMap();

    /** Returns array of elements who share given dof */
    Array<int> NeighborsWithDof(int e, int dof);
};

void CreateReferenceMesh(Mesh &reference_mesh,int dim);