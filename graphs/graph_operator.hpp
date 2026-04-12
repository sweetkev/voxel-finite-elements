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

    // Returns the GraphOperator corresponding the the next coarsened level of the graph/fes.
    GraphOperator Coarsen();

    SparseMatrix &GetMatrix() { return A; }

    /** Builds prolongartion matrix from coarse space to fine space using the 
    *   local prolongation defined on the reference element.
    */
    SparseMatrix CreateProlongation(DenseTensor &local_prolongation,
                                    Graph &fine_graph,
                                    const std::vector<int> &fine_broken_to_true_dof,
                                    const std::vector<std::set<int>> fine_dof_groups);

    const std::vector<int> &GetBrokenToTrueDofMap() const { return broken_to_true_dof; }

    const std::vector<std::set<int>> &GetDofGroups() const { return dof_groups; }

private:
    DenseMatrix A_ref;
    FiniteElementSpace reference_fes;

    // Matrix whose (i,j) entry denotes the local dof number on reference
    // element j which is identified with the local dof i on
    // the central reference element.
    // Is there a way to do this with a table instead of a DenseMatrix? Or
    // would a SparseMatrix be more appropriate?
    DenseMatrix local_to_neighbor_dof_map;

    // Maps each broken dof index to the true dof index.
    // Size is (ne * dofs_per_element).
    std::vector<int> broken_to_true_dof;

    SparseMatrix A;
    Graph graph;

    // Vector with size equal to number of true dofs whose entries contain the 
    // sets of broken dofs identified with the true dof.
    std::vector<std::set<int>> dof_groups;

    GraphOperator(DenseMatrix A_ref_, FiniteElementSpace &reference_fes_, 
                  Graph &graph_);

    /** Builds the matrix A_ref representing the stiffness matrix over the
     *  reference element.
     */
    void BuildARef();

    /** Creates the table whose ith row represents a local dof number
     *  and jth column represents a reference cell in the reference mesh
     */
    void CreateReferenceDofMap();

    /** Creates the vector with size equal to number of true dofs whose entries
     *  contain the sets of broken dofs identified with the true dof.
    */
    void CreateDofGroups();

    /** Returns the local dof number on neighbor element that corresponds
     *  to the local dof number e_dof on element e. If the dof is not shared
     *  between elements, returns -1
    */
    int GetNeighborDof(int e, int e_dof, int neighbor);

    /** Builds matrix A such that A = Lambda^T * hat{A} * Lambda. Here,
     *  hat{A} is the block diagonal matrix with blocks corresponding to the 
     *  local element matrix A_ref, and Lambda is the boolean matrix mapping
     *  broken dofs to the true dofs.
     */
    void BuildA(std::vector<std::set<int>> dof_groups);

    // Removes boundary dofs from the sparsematrix B
    void RemoveBoundaryDofs(SparseMatrix &B, 
                            Operator::DiagonalPolicy dpolicy);

};

/** Builds mesh which contains all neighbor information for a single element.
 *  For example in a 2D mesh, this is the 3x3 square. For a 3D mesh, it is
 *  The 3x3 cube.
 */
Mesh CreateReferenceMesh(int dim);
