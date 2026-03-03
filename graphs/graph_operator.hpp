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

    // Builds the coarse graph operator by coarsening the fine graph operator.
    GraphOperator Coarsen();

    SparseMatrix &GetMatrix() { return A; }

private:
    DenseMatrix A_ref;
    FiniteElementSpace reference_fes;

    // Matrix whose (i,j) entry denotes the local dof number on reference
    // element j which is identified with the local dof i on
    // the central reference element.
    DenseMatrix local_to_neighbor_dof_map;

    SparseMatrix A;
    Graph graph;

    GraphOperator(DenseMatrix &A_ref_, FiniteElementSpace &reference_fes_, 
                  Graph &graph_);

    /** Builds the matrix A_ref representing the stiffness matrix over the
     * reference element.
     */
    void BuildARef();

    /** Creates the table whose ith row represents a local dof number
     *  and jth column represents a reference cell in the reference mesh
     */
    void CreateDofMap();

    /** Creates the vector with size equal to number of true dofs whose entries
     *  contain the sets of broken dofs identified with the true dof.
    */
    std::vector<std::set<int>> CreateDofGroups();

    /** Returns the local dof number on neighbor element that corresponds
     *  to the local dof number e_dof on element e. If the dof is not shared
     *  between elements, returns -1
    */
    int GetNeighborDof(int e, int e_dof, int neighbor);

    /** Builds matrix A such that A = Lambda^T * hat{A} * Lambda. Here,
     * hat{A} is the block diagonal matrix with blocks corresponding to the 
     * local element matrix A_ref, and Lambda is the boolean matrix mapping the
     * broken dofs to the global dofs.
     */
    void BuildA(std::vector<std::set<int>> dof_groups);
};

/** Builds mesh which contains all neighbor information for a single element.
 *  For example in a 2D mesh, this is the 3x3 square. For a 3D mesh, it is
 *  The 3x3 cube.
 */
Mesh CreateReferenceMesh(int dim);
