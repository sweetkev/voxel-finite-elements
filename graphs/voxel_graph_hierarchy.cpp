#include "voxel_graph_hierarchy.hpp"
#include "mfem.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"

// Zeros the rows of P for the fine essential true dofs and the columns for the
// coarse essential true dofs. Collective, because EliminateCols is.
void EliminateProlongationEssentialDofs(HypreParMatrix &P,
                                        const Array<int> &coarse_ess_true_dofs,
                                        const Array<int> &fine_ess_true_dofs)
{
    P.EliminateRows(fine_ess_true_dofs);
    delete P.EliminateCols(coarse_ess_true_dofs);
}

VoxelGraphHierarchy::VoxelGraphHierarchy(unique_ptr<VoxelGraph> fine_graph,
                    Array<int> &fine_ess_dofs,
                    int nlevels,
                    const CellOwnership &ownership_,
                    FiniteElementSpace &reference_fes_,
                    real_t h_, ParBilinearForm &a, HypreParMatrix &fine_A_)
                    : fine_A(fine_A_), ownership(ownership_),
                      reference_fes(reference_fes_), h(h_)
{
    CreateLocalProlongation(reference_fes.GetMesh()->Dimension());

    graphs.resize(nlevels);
    ess_true_dofs.resize(nlevels);
    prolongations.resize(nlevels - 1);
    graph_operators.resize(nlevels);

    graphs[nlevels - 1] = std::move(fine_graph);
    ess_true_dofs[nlevels - 1] = fine_ess_dofs;
    graph_operators[nlevels - 1] = make_unique<VoxelGraphOperator>(reference_fes, *graphs[nlevels - 1], h, a);

    // Coarsen the graph level by level. The coarsening of each coarse grid
    // cell happens on its owner.
    for (int i = nlevels - 2; i >= 0; --i)
    {
        graphs[i] = make_unique<VoxelGraph>(graphs[i + 1]->CoarsenGraph(ownership, i));
    }

    // Build the prolongations and coarse operators.
    for (int i = nlevels - 2; i >= 0; --i)
    {
        prolongations[i] = CreateProlongation(*graphs[i], *graphs[i + 1],
                                              ess_true_dofs[i],
                                              ess_true_dofs[i + 1]);
        EliminateProlongationEssentialDofs(*prolongations[i], ess_true_dofs[i],
                                           ess_true_dofs[i + 1]);

        graph_operators[i] = make_unique<VoxelGraphOperator>(
            reference_fes, *graphs[i], *graph_operators[i + 1],
            local_prolongation);
        graph_operators[i]->GetParMatrix().EliminateBC(ess_true_dofs[i],
                                                       Operator::DIAG_ONE);
    }
}

void VoxelGraphHierarchy::CreateLocalProlongation(int dim)
{
    Mesh element;
    Geometry::Type geom;
    if (dim == 2)
    {
        element = Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL);
    }
    else if (dim == 3)
    {
        element = Mesh::MakeCartesian3D(1, 1, 1, Element::HEXAHEDRON);
    }
    else
    {
        MFEM_ABORT("Unsupported dimension.");
    }
    geom = element.GetTypicalElementGeometry();

    element.UniformRefinement();

    // Create local prolongation by refining a single element.
    const CoarseFineTransformations &rtrans = element.GetRefinementTransforms();
    const DenseTensor &pmats = rtrans.point_matrices[geom];
    const int nmat = pmats.SizeK();

    const FiniteElement *fe = reference_fes.GetFE(0);
    const int ldof = fe->GetDof();

    IsoparametricTransformation isotr;
    isotr.SetIdentityTransformation(geom);

    local_prolongation.SetSize(ldof, ldof, nmat);

    for (int i = 0; i < nmat; i++)
    {
        isotr.SetPointMat(pmats(i));
        fe->GetLocalInterpolation(isotr, local_prolongation(i));
    }
}

unique_ptr<HypreParMatrix> VoxelGraphHierarchy::CreateProlongation(
    const VoxelGraph &coarse_graph, const VoxelGraph &fine_graph,
    Array<int> &coarse_ess_true_dofs, const Array<int> &fine_ess_true_dofs)
{
    const vector<int> &fine_to_coarse_element_map = fine_graph.GetFineToCoarseElementMap();
    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int fine_nbdofs = fine_graph.Size() * dofs_per_elem;

    HypreParMatrix &fine_Lambda = fine_graph.GetBrokenToTrueDofMatrix();
    HypreParMatrix &coarse_Lambda = coarse_graph.GetBrokenToTrueDofMatrix();

    // Broken prolongation P-hat from coarse broken dofs to fine broken dofs. The
    // rows of local fine element e hold the local prolongation of e's
    // position in its coarse element, in the columns of the coarse element's
    // broken dofs.
    Array<int> I(fine_nbdofs + 1);
    Array<HYPRE_BigInt> J;
    J.Reserve(fine_nbdofs * dofs_per_elem);
    Array<real_t> values;
    values.Reserve(fine_nbdofs * dofs_per_elem);
    I[0] = 0;
    for (int fine_e = 0; fine_e < fine_graph.Size(); ++fine_e)
    {
        const int coarse_e = fine_to_coarse_element_map[fine_e];
        const DenseMatrix &block =
            local_prolongation(ChildIndex(fine_graph.GetElementCell(fine_e)));

        for (int rf = 0; rf < dofs_per_elem; ++rf)
        {
            const int row = fine_e * dofs_per_elem + rf;
            for (int rc = 0; rc < dofs_per_elem; ++rc)
            {
                if (block(rf, rc) == 0.0) { continue; }
                J.Append((HYPRE_BigInt)coarse_e * dofs_per_elem + rc);
                values.Append(block(rf, rc));
            }
            I[row + 1] = J.Size();
        }
    }
    const HypreParMatrix P_broken(fine_Lambda.GetComm(), fine_nbdofs,
                                  fine_Lambda.GetGlobalNumRows(),
                                  coarse_Lambda.GetGlobalNumRows(),
                                  I.GetData(), J.GetData(), values.GetData(),
                                  fine_Lambda.RowPart(), coarse_Lambda.RowPart());

    // P = D^{-1} * Lambda_f^T * P-hat * Lambda_c
    // Lambda_f^T sums the rows of P-hat over the broken copies of
    // each fine true dof. Every copy interpolates the same continuous coarse
    // function, so the copies' rows are equal and D^{-1}, with D the number of
    // copies of each fine true dof, turns the sum into that common row.
    const unique_ptr<HypreParMatrix> fine_Lambda_T(fine_Lambda.Transpose());
    const unique_ptr<HypreParMatrix> LtP(ParMult(fine_Lambda_T.get(), &P_broken, true));
    unique_ptr<HypreParMatrix> P(ParMult(LtP.get(), &coarse_Lambda, true));

    Vector ones(fine_Lambda.Height());
    ones = 1.0;
    Vector inverse_copies(fine_Lambda.Width());
    fine_Lambda.MultTranspose(ones, inverse_copies);
    for (real_t &copies : inverse_copies) { copies = 1.0 / copies; }
    P->ScaleRows(inverse_copies);

    // Identify the coarse essential true dofs: a coarse true dof is essential
    // if it has a nonzero contribution to any fine essential true dof.
    Vector fine_ess_marker(P->Height());
    fine_ess_marker = 0.0;
    for (int ess_dof : fine_ess_true_dofs) { fine_ess_marker[ess_dof] = 1.0; }
    Vector coarse_ess_weight(P->Width());
    P->AbsMultTranspose(fine_ess_marker, coarse_ess_weight);
    coarse_ess_true_dofs.DeleteAll();
    for (int i = 0; i < coarse_ess_weight.Size(); ++i)
    {
        if (coarse_ess_weight[i] > 1e-12) { coarse_ess_true_dofs.Append(i); }
    }

    return P;
}

Mesh CreateReferenceMesh(int dim)
{
    if (dim == 2)
    {
        return Mesh::MakeCartesian2D(3, 3, Element::QUADRILATERAL, false,
                                     1.0, 1.0, false);
    }
    else if (dim == 3)
    {
        return Mesh::MakeCartesian3D(3, 3, 3, Element::HEXAHEDRON,
                                     1.0, 1.0, 1.0, false);
    }
    else
    {
        MFEM_ABORT("Unsupported dimension.");
    }
}
