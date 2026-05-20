#include "voxel_graph_hierarchy.hpp"
#include "mfem.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"

VoxelGraphHierarchy::VoxelGraphHierarchy(unique_ptr<VoxelGraph> fine_graph, 
                    int nlevels, 
                    const FiniteElementSpace &reference_fes_)
                    : reference_fes(reference_fes_)
{
    CreateLocalProlongation(reference_fes.GetMesh()->Dimension());

    graphs.push_back(fine_graph);
    for (int i = 1; i < nlevels; ++i)
    {
        graphs.push_back(make_unique<VoxelGraph>(graphs[i-1]->CoarsenGraph()));
        // graph_operators.push_back(make_unique<VoxelGraphOperator>(*graphs[i]));
        prolongations.push_back(make_unique<SparseMatrix>(CreateProlongation(*graphs[i], *graphs[i-1])));
    }
}

void VoxelGraphHierarchy::CreateLocalProlongation(int dim)
{
    Mesh element;
    Geometry::Type geom;
    if (dim == 2)
    {
        element = Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL);
        geom = Geometry::SQUARE;
    }
    else if (dim == 3)
    {
        element = Mesh::MakeCartesian3D(1, 1, 1, Element::HEXAHEDRON);
        geom = Geometry::CUBE;
    }
    else
    {
        MFEM_ABORT("Unsupported dimension.");
    }

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