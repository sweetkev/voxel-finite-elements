# pragma once

#include "mfem.hpp"
#include "voxel_graph.hpp"

class VoxelGraphOperator
{
public:
    // Constructor for finest level using the bilinear form a. The finest 
    // parallel operator is MFEM's assembled matrix.
    VoxelGraphOperator(FiniteElementSpace &reference_fes_, VoxelGraph &graph_, real_t h_, BilinearForm &a);

    // Constructor for coarser levels which uses the finer operator. graph_ is
    // the coarsening of fine_operator's graph, and local_prolongation(k) is
    // the prolongation from a coarse element to its child with ChildIndex k.
    VoxelGraphOperator(FiniteElementSpace &reference_fes_, VoxelGraph &graph_,
                       VoxelGraphOperator &fine_operator,
                       const DenseTensor &local_prolongation);

    // Parallel matrix A_par = Λᵀ Â Λ in the graph's true dof numbering.
    HypreParMatrix &GetParMatrix() { return *A_par; }

    // Broken dofs to true dofs.
    HypreParMatrix &GetParLambda() { return graph.GetBrokenToTrueDofMatrix(); }

    // Block diagonal matrix of the element matrices of the local elements.
    SparseMatrix &GetAHat() { return A_hat; }

    VoxelGraph &GetGraph() { return graph; }

private:
    real_t h;

    DenseMatrix A_ref;
    FiniteElementSpace &reference_fes;

    SparseMatrix A_hat;

    Array<HYPRE_BigInt> broken_offsets;
    unique_ptr<HypreParMatrix> A_hat_par;
    unique_ptr<HypreParMatrix> A_par;

    VoxelGraph &graph;

    // Builds the matrix A_ref representing the stiffness matrix over the
    // reference element.
    void BuildARef(BilinearForm &a);

    // Builds hat{A}, the block diagonal matrix with a copy of A_ref for each
    // local element.
    void BuildAHat();

    // Builds hat{A} on a coarse level. The element matrix of coarse element E
    // is
    //     A_E = sum_{f} P_fᵀ A_f P_f,
    // where A_f is the element matrix of fine element f and P_f the local
    // prolongation for f's position in E. Each term is computed
    // on the rank holding f and sent to the rank holding E.
    void BuildAHat(VoxelGraphOperator &fine_operator,
                   const DenseTensor &local_prolongation);

    // Assembles A_par = Λᵀ Â Λ, with Λ the graph's broken to true dof matrix.
    void ParallelAssemble();
};
