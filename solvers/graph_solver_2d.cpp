#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include "graph.hpp"

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
    PixelMesh fine_mesh(pgm_file);

    //set up problem
    H1_FECollection fec(order, fine_mesh.GetMesh().Dimension());
    FiniteElementSpace fine_fes(&fine_mesh.GetMesh(), &fec);

    Array<int> fine_ess_dofs;
    fine_fes.GetBoundaryTrueDofs(fine_ess_dofs);

    BilinearForm fine_a(&fine_fes);
    fine_a.AddDomainIntegrator(new DiffusionIntegrator);
    fine_a.Assemble();

    GridFunction x(&fine_fes);
    x = 0.0;

    ConstantCoefficient one(1.0);
    LinearForm b(&fine_fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    SparseMatrix fine_A;
    Vector X, B;
    fine_a.FormLinearSystem(fine_ess_dofs, x, b, fine_A, X, B);

    //create multigrid operator
    Graph fine_graph(fine_fes);

    //solve problem

    fine_mesh.GetMesh().Save("mesh.mesh");

}