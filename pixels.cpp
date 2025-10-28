#include "ppm.hpp"
#include "mfem.hpp"
#include "pixel_mesh.hpp"
#include <fstream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

int main(int argc, char *argv[]) {
    //parse options
    std::string pgm_file = "pgm_files/small.pgm";

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&pgm_file, "-f", "--file", "pgm file to use");
    args.ParseCheck();

    //make mesh from file
    PixelMesh fine_mesh(pgm_file);

    //coarsen mesh
    PixelMesh coarse_mesh = fine_mesh.CoarsenMesh();

    /*

    //Define finite elment space on fine mesh
    H1_FECollection fec(1,2);
    FiniteElementSpace fespace(&fine_mesh.GetMesh(), &fec);

    //Get boundary dofs to enforce boundary conditions
    Array<int> boundary_dofs;
    fespace.GetBoundaryTrueDofs(boundary_dofs);

    //setup gridfunction x
    GridFunction x(&fespace);
    x = 0.0;

    //setup linear form b
    ConstantCoefficient one(1.0);
    LinearForm b(&fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    //setup multigrid
    Multigrid m;

    //solve
    SparseMatrix A;
    Vector B,X;
    m.FormLinearSystem(boundary_dofs, x, b, A, X, B);
    m.RecoverFEMSolution(X, b, x);

    */
}