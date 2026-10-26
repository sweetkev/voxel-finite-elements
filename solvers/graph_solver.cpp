#include "mfem.hpp"
#include "voxel_image.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"
#include "voxel_graph_hierarchy.hpp"

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
    // Initialize MPI and hypre.
    Mpi::Init();
    Hypre::Init();

    // Parse options
    string file = "../pgm_files/australia.pgm";
    int nlevels = -1;
    int order = 1;
    bool dirichlet = true;
    string mesh_file = "";

    OptionsParser args(argc, argv);
    args.AddOption(&file, "-f", "--file", "image file to use");
    args.AddOption(&mesh_file, "-m", "--mesh", "mesh file to use instead of image");
    args.AddOption(&nlevels, "-nl", "--nlevels", "number of multigrid levels");
    args.AddOption(&order, "-o", "--order", "polynomial order");
    args.AddOption(&dirichlet, "-d", "--dirichlet", "-n", "--neumann",
                   "Dirichlet or Neumann BCs");
    args.ParseCheck();

    // Generate fine mesh from image or mesh file
    unique_ptr<VoxelImage> image_ptr;
    if (mesh_file.empty())
    {
        image_ptr = make_unique<VoxelImage>(file);
    }
    else
    {
        image_ptr = make_unique<VoxelImage>(Mesh(mesh_file));
    }
    VoxelImage &image = *image_ptr;
    Mesh mesh = image.MakeMesh();

    const real_t h = mesh.GetElementSize(0, 0);

    // If nlevels was not specified, coarsen until smallest dimension is 1
    if (nlevels < 0)
    {
        Vector bb_min;
        Vector bb_max;
        mesh.GetBoundingBox(bb_min, bb_max, 0);
        Vector lengths(bb_max);
        lengths -= bb_min;
        real_t n = lengths.Min() / h;

        nlevels = 0;
        while (pow(2, nlevels) < n) { ++nlevels; }
        // One more level than coarsenings
        nlevels += 1;
        
        if(Mpi::Root()) { cout << "new nlevels: " << nlevels << "\n"; }
    }

    // Set up problem in parallel
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    H1_FECollection fec(order, mesh.Dimension());
    // FiniteElementSpace fes(&mesh, &fec);
    ParFiniteElementSpace pfes(&pmesh, &fec);

    Array<int> ess_dofs;
    if (dirichlet)
    {
        pfes.GetBoundaryTrueDofs(ess_dofs);
    }

    ParBilinearForm a(&pfes);
    a.AddDomainIntegrator(new DiffusionIntegrator);
    a.AddDomainIntegrator(new MassIntegrator);
    a.Assemble();

    ParGridFunction x(&pfes);
    x = 0.0;

    ConstantCoefficient one(1.0);
    ParLinearForm b(&pfes);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    HypreParMatrix A;
    Vector X, B;
    a.FormLinearSystem(ess_dofs, x, b, A, X, B);

    //create reference mesh for fespace
    Mesh reference_mesh = CreateReferenceMesh(mesh.Dimension());
    H1_FECollection reference_fec(order, reference_mesh.Dimension());
    FiniteElementSpace reference_fes(&reference_mesh, &reference_fec);

    // Create multigrid hierarchy
    VoxelGraphHierarchy graph_hierarchy(
        make_unique<VoxelGraph>(pfes, image, reference_fes), ess_dofs,
        nlevels, reference_fes, h, a);

    Array<Operator*> operators(nlevels);
    Array<Solver*> smoothers(nlevels);
    Array<Operator*> prolongations(nlevels - 1);
    Array<bool> own_operators(nlevels);
    Array<bool> own_smoothers(nlevels);
    Array<bool> own_prolongations(nlevels - 1);

    own_operators = false;
    own_smoothers = true;
    own_prolongations = false;

    for (int level = 0; level < nlevels - 1; ++level)
    {
        prolongations[level] = graph_hierarchy.GetProlongations()[level].get();
    }

    // TODO: Extend to n ranks
    // Wrap each level's SparseMatrix in a HypreParMatrix (single rank only:
    // the whole matrix is the diagonal block).
    vector<array<HYPRE_BigInt, 2>> row_starts(nlevels);
    vector<unique_ptr<HypreParMatrix>> par_operators(nlevels);

    // Set operators and smoothers for each level
    for (int level = 0; level < nlevels; level++)
    {
        // Convert SparseMatrix from graph code to ParMatrix
        SparseMatrix &current_A = graph_hierarchy.GetGraphOperators()[level]->GetMatrix();
        row_starts[level] = {0, current_A.Height()};
        par_operators[level] = make_unique<HypreParMatrix>(
            MPI_COMM_WORLD, current_A.Height(), row_starts[level].data(), &current_A);

        HypreParMatrix *A = par_operators[level].get();
        HYPRE_BigInt ndofs = A->GetGlobalNumRows();
        if (Mpi::Root()) { std::cout << "Level " << level << ": " << ndofs << " DOFs\n"; }
        operators[level] = A;

        if (level == 0)
        {
            // Coarsest level
            HypreBoomerAMG *amg = new HypreBoomerAMG(*A);
            amg->SetPrintLevel(0);
            smoothers[level] = amg;
        }
        else
        {
            // Finer levels: smoother
            HypreSmoother *s = new HypreSmoother(*A, HypreSmoother::l1GS);
            s->SetOperatorSymmetry(true);
            smoothers[level] = s;
        }
    }


    Multigrid mg(operators, smoothers, prolongations, own_operators,
                 own_smoothers, own_prolongations);

    //solve problem
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    //cg.SetOperator(graph_hierarchy.GetGraphOperators()[nlevels-1]->GetMatrix());
    cg.SetOperator(*par_operators[nlevels-1]);
    cg.SetPreconditioner(mg);
    cg.Mult(B, X);

    a.RecoverFEMSolution(X,b,x);
    x.SaveAsOne("sol.gf");
    pmesh.SaveAsOne("mesh.mesh");
}
