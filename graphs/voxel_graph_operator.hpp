# pragma once

#include "mfem.hpp"
#include "voxel_graph.hpp"

class VoxelGraphOperator
{
public:
    // Constructor for finest level using the bilinear form a.
    VoxelGraphOperator(FiniteElementSpace &reference_fes_, VoxelGraph &graph_, real_t h_, BilinearForm &a);

    // Constructor for coarser levels which uses the finer operator.
    VoxelGraphOperator(FiniteElementSpace &reference_fes_, VoxelGraph &graph_, VoxelGraphOperator &fine_operator, SparseMatrix &P);

    SparseMatrix &GetMatrix() { return A; }

    SparseMatrix &GetAHat() { return A_hat; }

private:
    real_t h;

    DenseMatrix A_ref;
    FiniteElementSpace &reference_fes;

    SparseMatrix A;
    SparseMatrix A_hat;

    VoxelGraph &graph;

    // Builds the matrix A_ref representing the stiffness matrix over the
    // reference element.
    void BuildARef(BilinearForm &a);

    // Builds matrix A such that A = Lambda^T * hat{A} * Lambda. Here,
    // hat{A} is the block diagonal matrix with blocks corresponding to the
    // local element matrix A_ref, and Lambda is the boolean matrix mapping
    // broken dofs to the true dofs.
    void BuildA(const std::vector<std::set<int>> dof_groups, 
                const std::vector<int> &broken_to_true_dof);

    // Builds the matrix A such that A = Lambda^T * hat{A} * Lambda. Here,
    // hat{A} is the block diagonal matrix with blocks corresponding to the
    // local element matrix A_e constructed from the hat{A} matrix from the 
    // finer level via 
    // A_e = (P_e)^T * (fine_Lambda_e)^T * fine_A_hat * fine_Lambda_e * P_e 
    // Lambda is the boolean matrix mapping broken dofs to the true dofs. 
    void BuildA(const std::vector<std::set<int>> dof_groups, 
                const std::vector<int> &broken_to_true_dof,
                VoxelGraphOperator &fine_operator,
                SparseMatrix &P);
};