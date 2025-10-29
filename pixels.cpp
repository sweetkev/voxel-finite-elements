#include "ppm.hpp"
#include "mfem.hpp"
#include "pixel_mesh.hpp"
#include <fstream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

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

   const int order = 1;

   H1_FECollection fec(order, fine_mesh.GetMesh().Dimension());
   FiniteElementSpace fine_fes(&fine_mesh.GetMesh(), &fec);
   FiniteElementSpace coarse_fes(&coarse_mesh.GetMesh(), &fec);

   const int nlevels = 2;

   Array<int> fine_ess_dofs;
   fine_fes.GetBoundaryTrueDofs(fine_ess_dofs);

   BilinearForm fine_a(&fespace);
   fine_a.AddDomainIntegrator(new DiffusionIntegrator);
   fine_a.Assemble();
   SparseMatrix fine_A;
   fine_a.FormSystemMatrix(fine_ess_dofs, fine_A);

   Array<int> coarse_ess_dofs;
   coarse_fes.GetBoundaryTrueDofs(coarse_ess_dofs);

   BilinearForm coarse_a(&fespace);
   coarse_a.AddDomainIntegrator(new DiffusionIntegrator);
   coarse_a.Assemble();
   SparseMatrix coarse_A;
   coarse_a.FormSystemMatrix(coarse_ess_dofs, coarse_A);

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

   Multigrid mg(operators, smoothers, prolongations, own_operators, own_smoothers,
                own_prolongations);

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
