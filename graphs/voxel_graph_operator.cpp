#include "mfem.hpp"
#include "voxel_graph_operator.hpp"
#include "voxel_graph.hpp"

VoxelGraphOperator::VoxelGraphOperator(FiniteElementSpace &reference_fes_,
                                         VoxelGraph &graph_, real_t h_, BilinearForm &a)
    : h(h_), reference_fes(reference_fes_), graph(graph_)
{
    // Build A_ref
    BuildARef(a);

    // Build A matrix
    BuildA(graph.NumLocalDofs(), graph.GetBrokenToLocalDofMap());
}

VoxelGraphOperator::VoxelGraphOperator(FiniteElementSpace &reference_fes_,
                                       VoxelGraph &graph_,
                                       VoxelGraphOperator &fine_operator,
                                       SparseMatrix &P)
                                       : reference_fes(reference_fes_), graph(graph_)
{
    BuildA(graph.NumLocalDofs(), graph.GetBrokenToLocalDofMap(), fine_operator, P);
}

void VoxelGraphOperator::ParallelAssemble(ParFiniteElementSpace &pfes)
{
    MPI_Comm comm = pfes.GetComm();

    // Row partition of the broken dofs (one block of rows per rank).
    HYPRE_BigInt local_nbdofs = A_hat.Height();
    HYPRE_BigInt global_nbdofs = 0;
    MPI_Allreduce(&local_nbdofs, &global_nbdofs, 1, HYPRE_MPI_BIG_INT,
                  MPI_SUM, comm);
    Array<HYPRE_BigInt> *offsets[1] = { &broken_offsets };
    pfes.GetParMesh()->GenerateOffsets(1, &local_nbdofs, offsets);

    // Â is block diagonal over elements, so it has no off-rank couplings.
    A_hat_par = make_unique<HypreParMatrix>(comm, global_nbdofs,
                                            broken_offsets.GetData(), &A_hat);

    // Λ_f = (broken -> local dof) * (local dof -> true dof)
    Lambda_local_par = make_unique<HypreParMatrix>(
        comm, global_nbdofs, pfes.GlobalVSize(), broken_offsets.GetData(),
        pfes.GetDofOffsets(), &Lambda);
    Lambda_par.reset(ParMult(Lambda_local_par.get(), pfes.Dof_TrueDof_Matrix()));

    A_par.reset(RAP(A_hat_par.get(), Lambda_par.get()));
}

void VoxelGraphOperator::BuildARef(BilinearForm &a)
{
    Array<BilinearFormIntegrator*> domain_integs = *a.GetDBFI();

    IsoparametricTransformation T;
    T.SetIdentityTransformation(reference_fes.GetMesh()->GetTypicalElementGeometry());

    auto &PM = T.GetPointMat();
    PM *= h;

    A_ref.SetSize(reference_fes.GetFE(0)->GetDof(), reference_fes.GetFE(0)->GetDof());
    A_ref = 0.0;

    for (BilinearFormIntegrator *integ : domain_integs)
    {
        DenseMatrix A_temp;
        integ->AssembleElementMatrix(*reference_fes.GetTypicalFE(), T, A_temp);
        A_ref += A_temp;
    }
}

void VoxelGraphOperator::BuildA(int nldofs,
                                const std::vector<int> &broken_to_local_dof)
{
    // Build A using matrix. Lambda(i,j) = 1 if broken dof i is identified with
    // local dof j, and 0 otherwise.

    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int ne = graph.Size();
    const int nbdofs = dofs_per_elem * ne;

    Lambda = SparseMatrix(nbdofs, nldofs);
    for (int bdof = 0; bdof < nbdofs; ++bdof)
    {
        int dof = broken_to_local_dof[bdof];
        Lambda.Set(bdof, dof, 1.0);
    }

    // Build hat{A} matrix.
    A_hat = SparseMatrix(nbdofs, nbdofs);

    for (int e = 0; e < ne; ++e)
    {
        Array<int> block_indices(dofs_per_elem);
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            block_indices[i] = e * dofs_per_elem + i;
        }
        A_hat.AddSubMatrix(block_indices, block_indices, A_ref);
    }

    // Finalize matrices before RAP operation
    Lambda.Finalize();
    A_hat.Finalize();

    // Compute A = Lambda^T * hat{A} * Lambda
    std::unique_ptr<SparseMatrix> A_ptr(RAP(Lambda, A_hat, Lambda));
    A.Swap(*A_ptr);
}

void VoxelGraphOperator::BuildA(int nldofs,
                                const std::vector<int> &broken_to_local_dof,
                                VoxelGraphOperator &fine_operator,
                                SparseMatrix &P)
{
    // Build A using matrix. Lambda(i,j) = 1 if broken dof i is identified with
    // local dof j, and 0 otherwise.

    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int ne = graph.Size();
    const int nbdofs = dofs_per_elem * ne;

    Lambda = SparseMatrix(nbdofs, nldofs);
    for (int bdof = 0; bdof < nbdofs; ++bdof)
    {
        int dof = broken_to_local_dof[bdof];
        Lambda.Set(bdof, dof, 1.0);
    }

    // Build hat{A} matrix.
    A_hat = SparseMatrix(nbdofs, nbdofs);

    const std::vector<int> &fine_broken_to_local_dof = 
        fine_operator.GetGraph().GetBrokenToLocalDofMap();

    for (int e = 0; e < ne; ++e)
    {
        Array<int> coarse_indices(dofs_per_elem);
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            coarse_indices[i] = e * dofs_per_elem + i;
        }

        const std::set<int> &contained_fine_elements = fine_operator.GetGraph().GetCoarseToFineElementsMap()[e];
        Array<int> fine_broken_dofs;
        for (int fe : contained_fine_elements)
        {
            for (int i = 0; i < dofs_per_elem; ++i)
            {
                fine_broken_dofs.Append(fe * dofs_per_elem + i);
            }
        }

        Array<int> fine_local_dofs;
        for (int i = 0; i < fine_broken_dofs.Size(); ++i)
        {
            int fine_local_dof = fine_broken_to_local_dof[fine_broken_dofs[i]];
            bool found = false;
            for (int dof : fine_local_dofs)
            {
                if (dof == fine_local_dof) { found = true; break; }
            }
            if (!found) { fine_local_dofs.Append(fine_local_dof); }
        }

        Array<int> coarse_local_dofs;
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            int coarse_broken_dof = e * dofs_per_elem + i;
            int coarse_local_dof = broken_to_local_dof[coarse_broken_dof];
            bool found = false;
            for (int dof : coarse_local_dofs)
            {
                if(dof == coarse_local_dof) { found = true; break; }
            }
            if (!found) { coarse_local_dofs.Append(coarse_local_dof); }
        }

        DenseMatrix local_A_hat(fine_broken_dofs.Size(), fine_broken_dofs.Size());
        fine_operator.GetAHat().GetSubMatrix(fine_broken_dofs, fine_broken_dofs, local_A_hat);

        DenseMatrix local_Lambda(fine_broken_dofs.Size(), fine_local_dofs.Size());
        fine_operator.GetLambda().GetSubMatrix(fine_broken_dofs, fine_local_dofs, local_Lambda);

        DenseMatrix local_P(fine_local_dofs.Size(), coarse_local_dofs.Size());
        P.GetSubMatrix(fine_local_dofs, coarse_local_dofs, local_P);

        DenseMatrix temp;
        DenseMatrix A_e;
        RAP(local_A_hat, local_Lambda, temp);
        RAP(temp, local_P, A_e);

        A_hat.AddSubMatrix(coarse_indices, coarse_indices, A_e);
    }

    // Finalize matrices before RAP operation
    Lambda.Finalize();
    A_hat.Finalize();

    // Compute A = Lambda^T * hat{A} * Lambda
    std::unique_ptr<SparseMatrix> A_ptr(RAP(Lambda, A_hat, Lambda));
    A.Swap(*A_ptr);
}

