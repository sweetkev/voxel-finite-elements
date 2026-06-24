#include "voxel_graph_hierarchy.hpp"
#include "mfem.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"

void EliminateProlongationEssentialDOFs(SparseMatrix &P,
                                        const Array<int> &coarse_ess_dofs,
                                        const Array<int> &fine_ess_dofs)
{
    for (const int f : fine_ess_dofs)
    {
        P.EliminateRow(f, Operator::DIAG_ZERO);
    }
    Array<int> coarse_ess_marker(P.Width());
    coarse_ess_marker = 0;
    for (int c : coarse_ess_dofs)
    {
        coarse_ess_marker[c] = 1;
    }
    P.EliminateCols(coarse_ess_marker);
}

VoxelGraphHierarchy::VoxelGraphHierarchy(unique_ptr<VoxelGraph> fine_graph,
                    Array<int> &fine_ess_dofs,
                    int nlevels,
                    FiniteElementSpace &reference_fes_,
                    real_t h_, BilinearForm &a)
                    : reference_fes(reference_fes_), h(h_)
{
    CreateLocalProlongation(reference_fes.GetMesh()->Dimension());

    graphs.resize(nlevels);
    ess_dofs.resize(nlevels);
    prolongations.resize(nlevels - 1);
    graph_operators.resize(nlevels);

    graphs[nlevels - 1] = std::move(fine_graph);
    ess_dofs[nlevels - 1] = fine_ess_dofs;
    graph_operators[nlevels - 1] = make_unique<VoxelGraphOperator>(reference_fes, *graphs[nlevels - 1], h, a);
    // Eliminate BCs from the finest level operator.
    for (const int ess_dof : fine_ess_dofs)
    {
        graph_operators[nlevels - 1]->GetMatrix().EliminateRowCol(ess_dof, Operator::DIAG_ONE);
    }

    for (int i = nlevels - 2; i >= 0; --i)
    {
        graphs[i] = make_unique<VoxelGraph>(graphs[i + 1]->CoarsenGraph());
        prolongations[i] = make_unique<SparseMatrix>(CreateProlongation(*graphs[i], *graphs[i+1], ess_dofs[i], ess_dofs[i+1]));
        graph_operators[i] = make_unique<VoxelGraphOperator>(reference_fes, *graphs[i], *graph_operators[i+1], *prolongations[i]);

        // Eliminate BCs from the operator and the prolongation matrix.
        for (const int ess_dof : ess_dofs[i])
        {
            graph_operators[i]->GetMatrix().EliminateRowCol(ess_dof, Operator::DIAG_ONE);
        }
        EliminateProlongationEssentialDOFs(*prolongations[i], ess_dofs[i], ess_dofs[i+1]);
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

SparseMatrix VoxelGraphHierarchy::CreateProlongation(const VoxelGraph &coarse_graph,
                                        const VoxelGraph &fine_graph,
                                        Array<int> &coarse_ess_dofs,
                                        const Array<int> &fine_ess_dofs)
{
    static constexpr int map_from_lex[4] = {0, 1, 3, 2};

    const vector<set<int>> &coarse_dof_groups = coarse_graph.GetDofGroups();
    const vector<int> &coarse_broken_to_true_dof = coarse_graph.GetBrokenToTrueDofMap();
    const vector<set<int>> &fine_dof_groups = fine_graph.GetDofGroups();
    const vector<int> &fine_broken_to_true_dof = fine_graph.GetBrokenToTrueDofMap();
    const vector<int> &fine_to_coarse_element_map = fine_graph.GetFineToCoarseElementMap();

    const int dofs_per_elem = reference_fes.GetFE(0)->GetDof();
    const int num_coarse_dofs = coarse_dof_groups.size();
    const int num_fine_dofs = fine_dof_groups.size();

    SparseMatrix P(num_fine_dofs, num_coarse_dofs);

    for (int fine_e = 0; fine_e < fine_graph.Size(); ++fine_e)
    {
        const int coarse_e = fine_to_coarse_element_map[fine_e];
        if (coarse_e < 0) { continue; }

        Coord coarse_coord = coarse_graph.GetElementCoord(coarse_e);
        Coord fine_coord = fine_graph.GetElementCoord(fine_e);

        // TODO: Generalize to n-dimensions. Currently assumes 2D.
        const int dx = fine_coord[0] - 2 * coarse_coord[0];
        const int dy = fine_coord[1] - 2 * coarse_coord[1];
        if (dx < 0 || dx > 1 || dy < 0 || dy > 1) { MFEM_ABORT(""); }

        const int local_prolongation_index = map_from_lex[dx + 2 * dy];

        Array<int> rows(dofs_per_elem);
        Array<int> cols(dofs_per_elem);
        for (int rf = 0; rf < dofs_per_elem; ++rf)
        {
            const int fine_tdof = fine_broken_to_true_dof[fine_e * dofs_per_elem + rf];
            rows[rf] = fine_tdof;
        }
        for (int rc = 0; rc < dofs_per_elem; ++rc)
        {
            const int coarse_tdof = coarse_broken_to_true_dof[coarse_e * dofs_per_elem + rc];
            cols[rc] = coarse_tdof;
        }
        P.SetSubMatrix(rows, cols, local_prolongation(local_prolongation_index));
    }
    P.Finalize();

    // Identify the coarse essential DOFs: a DOF is a coarse essential DOF if it
    // has a nonzero contribution to any fine essential DOF.
    coarse_ess_dofs.DeleteAll();
    if(!fine_ess_dofs.IsEmpty())
    {
        vector<bool> coarse_ess_dof_added(num_coarse_dofs, false);
        for (int ess_dof : fine_ess_dofs)
        {
            Array<int> cols;
            Vector row;
            P.GetRow(ess_dof, cols, row);
            for (int i = 0; i < cols.Size(); ++i)
            {
                const int col = cols[i];
                if(abs(row[i]) > 1e-12 && !coarse_ess_dof_added[col])
                {
                    coarse_ess_dofs.Append(col);
                    coarse_ess_dof_added[col] = true;
                }
            }
        }
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
