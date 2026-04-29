#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include "graph.hpp"
#include "graph_operator.hpp"

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{

    //parse options
    string pgm_file = "../pgm_files/australia.pgm";
    int nlevels = -1;
    int order = 1;

    OptionsParser args(argc, argv);
    args.AddOption(&pgm_file, "-f", "--file", "pgm file to use");
    args.AddOption(&nlevels, "-nl", "--nlevels", "number of multigrid levels");
    args.AddOption(&order, "-o", "--order", "polynomial order");
    args.ParseCheck();

    //generate fine mesh from pgm file
    PixelImage image(pgm_file);
    PixelMesh mesh(image);

    //set up problem
    H1_FECollection fec(order, mesh.Dimension());
    FiniteElementSpace fes(&mesh.GetMesh(), &fec);

    Array<int> ess_dofs;
    // fes.GetBoundaryTrueDofs(ess_dofs);

    BilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator);
    a.AddDomainIntegrator(new MassIntegrator);
    a.Assemble();

    GridFunction x(&fes);
    x = 0.0;

    ConstantCoefficient one(1.0);
    LinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    SparseMatrix A;
    Vector X, B;
    a.FormLinearSystem(ess_dofs, x, b, A, X, B);

    //create multigrid operator
    //create reference mesh
    Mesh reference_mesh = CreateReferenceMesh(mesh.GetMesh().Dimension());
    H1_FECollection reference_fec(order, reference_mesh.Dimension());
    FiniteElementSpace reference_fes(&reference_mesh, &reference_fec);

    // Create fine graph
    Graph fine_graph(fes, image);

    const real_t h = mesh.GetMesh().GetElementSize(0, 0);

    // Create a GraphOperator for the fine graph to get the broken->true dof map
    GraphOperator fine_graph_operator(reference_fes, fine_graph, h);

    // Create coarse graph and its corresponding operator
    GraphOperator coarse_graph_operator = fine_graph_operator.Coarsen();
    SparseMatrix &fine_A = fine_graph_operator.GetMatrix();
    SparseMatrix &coarse_A = coarse_graph_operator.GetMatrix();

    {
        std::ofstream f("A.txt");
        A.PrintMatlab(f);
    }
    {
        std::ofstream f("coarse_A.txt");
        coarse_A->PrintMatlab(f);
    }

    // Create prolongation operator from coarse to fine space
    // Ensure reference element fespace matches the others used
    Mesh reference_element = Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL);
    H1_FECollection single_fec(order, reference_element.Dimension());
    FiniteElementSpace single_fes(&reference_element, &single_fec);

    reference_element.UniformRefinement();
    single_fes.Update();

    const Geometry::Type geom = Geometry::SQUARE;

    const CoarseFineTransformations &rtrans = reference_element.GetRefinementTransforms();
    const DenseTensor &pmats = rtrans.point_matrices[geom];
    const int nmat = pmats.SizeK();

    const FiniteElement *fe = single_fes.GetFE(0);
    const int ldof = fe->GetDof();

    IsoparametricTransformation isotr;
    isotr.SetIdentityTransformation(geom);

    DenseTensor local_prolongation(ldof, ldof, nmat);

    for (int i = 0; i < nmat; i++)
    {
        isotr.SetPointMat(pmats(i));
        fe->GetLocalInterpolation(isotr, local_prolongation(i));
    }

    SparseMatrix *P = new SparseMatrix(coarse_graph_operator.CreateProlongation(
        local_prolongation,
        fine_graph,
        fine_graph_operator.GetBrokenToTrueDofMap(),
        fine_graph_operator.GetDofGroups(),
        coarse_graph_operator.GetGraph().GetGraphLabeling()));

    {
        std::ofstream f("P.txt");
        P->PrintMatlab(f);
    }

    // Create multigrid hierarchy
    Array<Operator*> operators(2);
    Array<Solver*> smoothers(2);
    Array<Operator*> prolongations(1);
    Array<bool> own_operators(2);
    Array<bool> own_smoothers(2);
    Array<bool> own_prolongations(1);

    operators[0] = &coarse_A;
    operators[1] = &fine_A;

    smoothers[0] = new UMFPackSolver(coarse_A);
    smoothers[1] = new GSSmoother(fine_A);

    prolongations[0] = P;

    own_operators[0] = false;
    own_operators[1] = false;
    own_smoothers = true;
    own_prolongations = true;


    Multigrid mg(operators, smoothers, prolongations, own_operators,
                 own_smoothers, own_prolongations);

    //solve problem
    CGSolver cg;
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    cg.SetOperator(fine_A);
    cg.SetPreconditioner(mg);
    cg.Mult(B, X);

    a.RecoverFEMSolution(X,b,x);
    x.Save("sol.gf");
    mesh.GetMesh().Save("mesh.mesh");
}
