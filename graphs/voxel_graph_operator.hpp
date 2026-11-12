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

    // For the finest level, assembles A_par = Λ_fᵀ Â Λ_f in pfes's true dof
    // numbering, where Λ_f = Λ * pfes.Dof_TrueDof_Matrix() maps broken dofs
    // to global true dofs.
    void ParallelAssemble(ParFiniteElementSpace &pfes);

    // Local matrix A = Λᵀ Â Λ. On the finest level its rows/columns are
    // pfes's local dofs, which coincide with true dofs only on one rank.
    SparseMatrix &GetMatrix() { return A; }

    HypreParMatrix &GetParMatrix()
    { return *A_par; }

    HypreParMatrix &GetParLambda()
    { return *Lambda_par; }

    SparseMatrix &GetAHat() { return A_hat; }

    SparseMatrix &GetLambda() { return Lambda; }

    VoxelGraph &GetGraph() { return graph; }

private:
    real_t h;

    DenseMatrix A_ref;
    FiniteElementSpace &reference_fes;

    SparseMatrix A;
    SparseMatrix A_hat;
    SparseMatrix Lambda;

    // Parallel versions (see ParallelAssemble). The block-diagonal wrappers
    // alias A_hat and Lambda, so those must outlive them.
    Array<HYPRE_BigInt> broken_offsets;
    unique_ptr<HypreParMatrix> A_hat_par;
    unique_ptr<HypreParMatrix> Lambda_local_par; // broken -> local dofs
    unique_ptr<HypreParMatrix> Lambda_par;       // broken -> true dofs
    unique_ptr<HypreParMatrix> A_par;

    VoxelGraph &graph;

    // Builds the matrix A_ref representing the stiffness matrix over the
    // reference element.
    void BuildARef(BilinearForm &a);

    // Builds matrix A such that A = Lambda^T * hat{A} * Lambda. Here,
    // hat{A} is the block diagonal matrix with blocks corresponding to the
    // local element matrix A_ref, and Lambda is the boolean matrix mapping
    // broken dofs to the true dofs.
    void BuildA(int ntdofs, const std::vector<int> &broken_to_true_dof);

    // Builds the matrix A such that A = Lambda^T * hat{A} * Lambda. Here,
    // hat{A} is the block diagonal matrix with blocks corresponding to the
    // local element matrix A_e constructed from the hat{A} matrix from the 
    // finer level via 
    // A_e = (P_e)^T * (fine_Lambda_e)^T * fine_A_hat * fine_Lambda_e * P_e 
    // Lambda is the boolean matrix mapping broken dofs to the true dofs. 
    void BuildA(int ntdofs, const std::vector<int> &broken_to_true_dof,
                VoxelGraphOperator &fine_operator,
                SparseMatrix &P);
};