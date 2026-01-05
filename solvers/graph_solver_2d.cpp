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
    PixelImage image(pgm_file);
    PixelMesh mesh(image);

    //set up problem
    H1_FECollection fec(order, mesh.GetMesh().Dimension());
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
    Graph fine_graph(fes, image);
    Graph coarse_graph = fine_graph.CoarsenGraph();

    //solve problem

    mesh.GetMesh().Save("mesh.mesh");

}