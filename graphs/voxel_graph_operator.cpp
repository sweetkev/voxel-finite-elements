#include "mfem.hpp"
#include "voxel_graph_operator.hpp"
#include "voxel_graph.hpp"
#include "cell_ownership.hpp"

VoxelGraphOperator::VoxelGraphOperator(FiniteElementSpace &reference_fes_,
                                         VoxelGraph &graph_, real_t h_, BilinearForm &a)
    : h(h_), reference_fes(reference_fes_), graph(graph_)
{
    BuildARef(a);
    BuildAHat();
    ParallelAssemble();
}

VoxelGraphOperator::VoxelGraphOperator(FiniteElementSpace &reference_fes_,
                                       VoxelGraph &graph_,
                                       VoxelGraphOperator &fine_operator,
                                       const DenseTensor &local_prolongation)
                                       : reference_fes(reference_fes_), graph(graph_)
{
    BuildAHat(fine_operator, local_prolongation);
    ParallelAssemble();
}

void VoxelGraphOperator::ParallelAssemble()
{
    HypreParMatrix &Lambda_par = graph.GetBrokenToTrueDofMatrix();
    MPI_Comm comm = Lambda_par.GetComm();
    int nranks;
    MPI_Comm_size(comm, &nranks);

    // Â has the row partition of the broken dofs, the rows of Λ.
    const int nparts = HYPRE_AssumedPartitionCheck() ? 2 : nranks + 1;
    broken_offsets.SetSize(nparts);
    for (int i = 0; i < nparts; ++i)
    {
        broken_offsets[i] = Lambda_par.RowPart()[i];
    }

    // Â is block diagonal over elements, so it has no off-rank couplings.
    A_hat_par = make_unique<HypreParMatrix>(comm, Lambda_par.GetGlobalNumRows(),
                                            broken_offsets.GetData(), &A_hat);

    A_par.reset(RAP(A_hat_par.get(), &Lambda_par));
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

void VoxelGraphOperator::BuildAHat()
{
    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int ne = graph.Size();
    const int nbdofs = dofs_per_elem * ne;

    A_hat = SparseMatrix(nbdofs, nbdofs);

    Array<int> block_indices(dofs_per_elem);
    for (int e = 0; e < ne; ++e)
    {
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            block_indices[i] = e * dofs_per_elem + i;
        }
        A_hat.AddSubMatrix(block_indices, block_indices, A_ref);
    }
    A_hat.Finalize();
}

void VoxelGraphOperator::BuildAHat(VoxelGraphOperator &fine_operator,
                                   const DenseTensor &local_prolongation)
{
    const VoxelGraph &fine_graph = fine_operator.GetGraph();
    const SparseMatrix &fine_A_hat = fine_operator.GetAHat();
    const vector<int> &fine_to_coarse = fine_graph.GetFineToCoarseElementMap();

    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int block_size = dofs_per_elem * dofs_per_elem;
    const int ne = graph.Size();
    const int nbdofs = dofs_per_elem * ne;

    MPI_Comm comm = graph.GetBrokenToTrueDofMatrix().GetComm();
    int nranks;
    MPI_Comm_size(comm, &nranks);

    // On the rank of each fine element f, compute P_fᵀ A_f P_f and send it,
    // with the global number of f's coarse element, to that element's rank.
    vector<vector<int>> send_elements(nranks);
    vector<vector<real_t>> send_blocks(nranks);
    Array<int> fine_dofs(dofs_per_elem);
    DenseMatrix A_f(dofs_per_elem);
    DenseMatrix C_f;
    for (int f = 0; f < fine_graph.Size(); ++f)
    {
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            fine_dofs[i] = f * dofs_per_elem + i;
        }
        fine_A_hat.GetSubMatrix(fine_dofs, fine_dofs, A_f);

        const DenseMatrix &P_f =
            local_prolongation(ChildIndex(fine_graph.GetElementCell(f)));
        RAP(A_f, P_f, C_f);

        const int coarse_global = fine_to_coarse[f];
        const int dest = graph.ElementRank(coarse_global);
        send_elements[dest].push_back(coarse_global);
        send_blocks[dest].insert(send_blocks[dest].end(), C_f.Data(),
                                 C_f.Data() + block_size);
    }

    vector<int> recv_counts;
    const vector<int> recv_elements = Exchange(comm, send_elements, recv_counts);
    vector<real_t> recv_blocks = Exchange(comm, send_blocks, recv_counts);

    // Sum the received terms into the coarse element matrices.
    A_hat = SparseMatrix(nbdofs, nbdofs);
    const int first_element = graph.GlobalElementNumber(0);
    Array<int> coarse_dofs(dofs_per_elem);
    for (size_t j = 0; j < recv_elements.size(); ++j)
    {
        const int e = recv_elements[j] - first_element;
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            coarse_dofs[i] = e * dofs_per_elem + i;
        }
        const DenseMatrix C(recv_blocks.data() + j * block_size,
                            dofs_per_elem, dofs_per_elem);
        A_hat.AddSubMatrix(coarse_dofs, coarse_dofs, C);
    }
    A_hat.Finalize();
}
