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
    BuildA(graph.GetDofGroups(), graph.GetBrokenToTrueDofMap());
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

void VoxelGraphOperator::BuildA(const std::vector<std::set<int>> dof_groups,
                                const std::vector<int> &broken_to_true_dof)
{
    // Build A using matrix. Lambda(i,j) = 1 if broken dof i is identified with
    // true dof j, and 0 otherwise.

    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int ne = graph.Size();
    const int nbdofs = dofs_per_elem * ne;
    const int ntdofs = dof_groups.size();

    SparseMatrix Lambda(nbdofs, ntdofs);
    for (int bdof = 0; bdof < nbdofs; ++bdof)
    {
        int tdof = broken_to_true_dof[bdof];
        Lambda.Set(bdof, tdof, 1.0);
    }

    // Build hat{A} matrix.
    SparseMatrix A_hat(nbdofs, nbdofs);

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

    // TODO: Eliminate boundary dofs from coarse operator
}
