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

    // If nlevels was not specified, coarsen until smallest dimension is 1
    if(nlevels < 0) {
      int width = mesh.GetWidth();
      int height = mesh.GetHeight();
      
      // Find the number of coarsenings until height or width is 1
      nlevels = max(ceil(log2(width)), ceil(log2(height)));
      cout << "new nlevels: " << nlevels << "\n";
    }

    //create reference mesh and compute local prolongation
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

    //create reference mesh for fespace
    Mesh reference_mesh = CreateReferenceMesh(mesh.GetMesh().Dimension());
    H1_FECollection reference_fec(order, reference_mesh.Dimension());
    FiniteElementSpace reference_fes(&reference_mesh, &reference_fec);

    const real_t h = mesh.GetMesh().GetElementSize(0, 0);

    // Create multigrid hierarchy using graph operators
    vector<unique_ptr<GraphOperator>> graph_operators(nlevels);

    // Create fine graph operator
    Graph fine_graph(fes, image);
    graph_operators[nlevels-1] = make_unique<GraphOperator>(reference_fes,
                                                              fine_graph,
                                                              h);

    // Create coarser levels
    for (int level = nlevels-2; level >= 0; --level)
    {
        graph_operators[level] = make_unique<GraphOperator>(
            graph_operators[level+1]->Coarsen());
    }

    // Create multigrid hierarchy
    Array<Operator*> operators(nlevels);
    Array<Solver*> smoothers(nlevels);
    Array<Operator*> prolongations(nlevels - 1);
    Array<bool> own_operators(nlevels);
    Array<bool> own_smoothers(nlevels);
    Array<bool> own_prolongations(nlevels - 1);

    own_operators = false;
    own_smoothers = true;
    own_prolongations = false;

    // Set operators and smoothers for each level
    for (int level = 0; level < nlevels; level++)
    {
        SparseMatrix *A = &graph_operators[level]->GetMatrix();
        operators[level] = A;
        
        if (level == 0)
        {
            // Coarsest level: direct solver
            smoothers[level] = new UMFPackSolver(*A);
        }
        else
        {
            // Finer levels: smoother
            smoothers[level] = new GSSmoother(*A);
        }
    }

    // Create prolongation operators between levels
    for (int level = 0; level < nlevels - 1; level++)
    {
        SparseMatrix *P = new SparseMatrix(
            graph_operators[level]->CreateProlongation(
                local_prolongation,
                graph_operators[level+1]->GetGraph(),
                graph_operators[level+1]->GetBrokenToTrueDofMap(),
                graph_operators[level+1]->GetDofGroups(),
                graph_operators[level]->GetGraph().GetGraphLabeling()));
        prolongations[level] = P;
        own_prolongations[level] = true;
    }


    Multigrid mg(operators, smoothers, prolongations, own_operators,
                 own_smoothers, own_prolongations);

    //solve problem
    CGSolver cg;
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    cg.SetOperator(graph_operators[nlevels-1]->GetMatrix());
    cg.SetPreconditioner(mg);
    cg.Mult(B, X);

    a.RecoverFEMSolution(X,b,x);
    x.Save("sol.gf");
    mesh.GetMesh().Save("mesh.mesh");
}
