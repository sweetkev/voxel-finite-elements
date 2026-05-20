#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"
#include "voxel_graph_hierarchy.hpp"

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

    //create reference mesh for fespace
    Mesh reference_mesh = CreateReferenceMesh(mesh.GetMesh().Dimension());
    H1_FECollection reference_fec(order, reference_mesh.Dimension());
    FiniteElementSpace reference_fes(&reference_mesh, &reference_fec);

    const real_t h = mesh.GetMesh().GetElementSize(0, 0);

    // Create multigrid hierarchy
    VoxelGraph fine_graph(fes, image, reference_fes);
    VoxelGraphHierarchy graph_hierarchy(make_unique<VoxelGraph>(fine_graph), nlevels, reference_fes, h);

    // Create multigrid hierarchy
    Array<Operator*> operators(nlevels);
    Array<Solver*> smoothers(nlevels);
    Array<Operator*> prolongations(nlevels - 1);
    Array<bool> own_operators(nlevels);
    Array<bool> own_smoothers(nlevels);
    Array<bool> own_prolongations(nlevels - 1);

    own_operators = false;
    own_smoothers = true;
    own_prolongations = true;

    prolongations = graph_hierarchy.GetProlongations();

    // Set operators and smoothers for each level
    for (int level = 0; level < nlevels; level++)
    {
        SparseMatrix *A = &graph_hierarchy.GetGraphOperators()[level]->GetMatrix();
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


    Multigrid mg(operators, smoothers, prolongations, own_operators,
                 own_smoothers, own_prolongations);

    //solve problem
    CGSolver cg;
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    cg.SetOperator(graph_hierarchy.GetGraphOperators()[nlevels-1]->GetMatrix());
    cg.SetPreconditioner(mg);
    cg.Mult(B, X);

    a.RecoverFEMSolution(X,b,x);
    x.Save("sol.gf");
    mesh.GetMesh().Save("mesh.mesh");
}
