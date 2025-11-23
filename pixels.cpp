#include "ppm.hpp"
#include "mfem.hpp"
#include "pixel_mesh.hpp"
#include <fstream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

void AddProlongationBCs(SparseMatrix &P, Array<int> &fine_ess_dofs);

int main(int argc, char *argv[])
{
   //parse options
   std::string pgm_file = "pgm_files/small.pgm";

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&pgm_file, "-f", "--file", "pgm file to use");
   args.ParseCheck();

   //make mesh from file
   PixelMesh fine_mesh(pgm_file);

   //coarsen mesh
   PixelMesh coarse_mesh = fine_mesh.CoarsenMesh();

   coarse_mesh.GetMesh().Save("coarse_mesh.mesh");

   const int order = 1;

   //setup finite element spaces
   H1_FECollection fec(order, fine_mesh.GetMesh().Dimension());
   FiniteElementSpace fine_fes(&fine_mesh.GetMesh(), &fec);
   FiniteElementSpace coarse_fes(&coarse_mesh.GetMesh(), &fec);

   Array<int> fine_ess_dofs;
   fine_fes.GetBoundaryTrueDofs(fine_ess_dofs);

   BilinearForm fine_a(&fine_fes);
   fine_a.AddDomainIntegrator(new DiffusionIntegrator);
   fine_a.Assemble();

   //initial guess
   GridFunction x(&fine_fes);
   x = 0.0;

   //right-hand side
   ConstantCoefficient one(1.0);
   LinearForm b(&fine_fes);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   SparseMatrix fine_A;
   Vector X, B;
   fine_a.FormLinearSystem(fine_ess_dofs, x, b, fine_A, X, B);

   Array<int> coarse_ess_dofs;
   coarse_fes.GetBoundaryTrueDofs(coarse_ess_dofs);

   BilinearForm coarse_a(&coarse_fes);
   coarse_a.AddDomainIntegrator(new DiffusionIntegrator);
   coarse_a.Assemble();
   SparseMatrix coarse_A;
   coarse_a.FormSystemMatrix(coarse_ess_dofs, coarse_A);

   //create multigrid operators
   const int nlevels = 2;
   Array<Operator*> operators(nlevels);
   Array<Solver*> smoothers(nlevels);
   Array<Operator*> prolongations(nlevels - 1);
   Array<bool> own_operators(nlevels);
   Array<bool> own_smoothers(nlevels);
   Array<bool> own_prolongations(nlevels - 1);

   own_operators = false;
   own_smoothers = false;
   own_prolongations = false;

   operators[0] = &coarse_A;
   operators[1] = &fine_A;

   UMFPackSolver coarse_solver(coarse_A);
   DSmoother fine_smoother(fine_A);
   smoothers[0] = &coarse_solver;
   smoothers[1] = &fine_smoother;

   SparseMatrix P = CreatePixelProlongation(coarse_mesh, coarse_fes, fine_mesh, fine_fes, fine_ess_dofs);
   prolongations[0] = &P;

   Multigrid mg(operators, smoothers, prolongations, own_operators, own_smoothers,
                own_prolongations);

   CGSolver cg;
   cg.SetRelTol(1e-12);
   cg.SetMaxIter(2000);
   cg.SetPrintLevel(1);
   cg.SetOperator(fine_A);
   cg.SetPreconditioner(mg);
   cg.Mult(B, X);

   fine_a.RecoverFEMSolution(X,b,x);
   x.Save("sol.gf");
   fine_mesh.GetMesh().Save("fine_mesh.mesh");
}
