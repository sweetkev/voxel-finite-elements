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
    fes.GetBoundaryTrueDofs(ess_dofs);

    BilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator);
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

    // Create coarse graph and its corresponding operator
    Graph coarse_graph = fine_graph.CoarsenGraph();
    
    GraphOperator coarse_graph_operator(reference_fes, coarse_graph);
    SparseMatrix coarse_A = coarse_graph_operator.GetMatrix();

    // Create prolongation operator from coarse to fine space
    // This is where I feel least comfortable. Below is the code that was
    // talked about several weeks ago.

    // Mesh mesh_ = Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL);

    // H1_FECollection fec(5, mesh_.Dimension());
    // FiniteElementSpace fespace(&mesh_, &fec);

    // mesh_.UniformRefinement();
    // fespace.Update();

    // const Geometry::Type geom = Geometry::SQUARE;

    // const CoarseFineTransformations &rtrans = mesh_.GetRefinementTransforms();
    // const DenseTensor &pmats = rtrans.point_matrices[geom];
    // const int nmat = pmats.SizeK();

    // const FiniteElement *fe = fespace.GetFE(0);
    // const int ldof = fe->GetDof();

    // IsoparametricTransformation isotr;
    // isotr.SetIdentityTransformation(geom);

    // DenseMatrix P(ldof, ldof);

    // for (int i = 0; i < nmat; i++)
    // {
    //     cout << "Matrix " << i << "\n";
    //     cout << "====================\n";

    //     isotr.SetPointMat(pmats(i));
    //     fe->GetLocalInterpolation(isotr, P);

    //     P.Print(cout, 4);

    //     cout << "\n\n";
    // }

    // // Create multigrid hierarchy
    // Array<Operator*> operators(2);
    // Array<Solver*> smoothers(2);
    // Array<Operator*> prolongations(1);
    // Array<bool> own_operators(2);
    // Array<bool> own_smoothers(2);
    // Array<bool> own_prolongations(1);

    // operators[0] = &coarse_A;
    // operators[1] = &A;

    // smoothers[0] = new UMFPackSolver(A);
    // smoothers[1] = new GSSmoother(coarse_A);

    // prolongations[0] = &P;

    // own_operators = false;
    // own_smoothers = false;
    // own_prolongations = false;

    // Multigrid mg(operators, smoothers, prolongations, own_operators, 
    //              own_smoothers, own_prolongations);

    // //solve problem
    // CGSolver cg;
    // cg.SetRelTol(1e-12);
    // cg.SetMaxIter(2000);
    // cg.SetPrintLevel(1);
    // cg.SetOperator(A);
    // cg.SetPreconditioner(mg);
    // cg.Mult(B, X);

    // a.RecoverFEMSolution(X,b,x);
    // x.Save("sol.gf");
    // mesh.GetMesh().Save("mesh.mesh");
}
