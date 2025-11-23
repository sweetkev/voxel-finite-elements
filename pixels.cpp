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
   std::string pgm_file = "pgm_files/australia.pgm";

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&pgm_file, "-f", "--file", "pgm file to use");
   args.ParseCheck();

   //make mesh from file
   PixelMesh fine_mesh(pgm_file);

   const int order = 1;

   H1_FECollection fec(order, fine_mesh.GetMesh().Dimension());
   FiniteElementSpace fine_fes(&fine_mesh.GetMesh(), &fec);

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

   Array<PixelMesh*> meshes(nlevels);
   Array<FiniteElementSpace*> fe_spaces(nlevels);
   Array<Array<int>*> ess_dofs(nlevels);

   meshes[nlevels-1] = &fine_mesh;
   fe_spaces[nlevels-1] = &fine_fes;
   ess_dofs[nlevels-1] = &fine_ess_dofs;

   operators[nlevels-1] = &fine_A;

   DSmoother fine_smoother(fine_A);
   smoothers[nlevels-1] = &fine_smoother;

   for(int level = nlevels-2; level >=0; --level)
   {
      PixelMesh coarse_mesh = meshes[level+1]->CoarsenMesh();
      meshes[level] = &coarse_mesh;

      FiniteElementSpace coarse_fes(&coarse_mesh.GetMesh(), &fec);
      fe_spaces[level] = &coarse_fes;

      Array<int> coarse_ess_dofs;
      coarse_fes.GetBoundaryTrueDofs(coarse_ess_dofs);
      ess_dofs[level] = &coarse_ess_dofs;

      BilinearForm coarse_a(&coarse_fes);
      coarse_a.AddDomainIntegrator(new DiffusionIntegrator);
      coarse_a.Assemble();
      SparseMatrix coarse_A;
      coarse_a.FormSystemMatrix(coarse_ess_dofs, coarse_A);

      operators[level] = &coarse_A;
      
      if(level == 0) {
         UMFPackSolver coarse_solver(coarse_A);
         smoothers[0] = &coarse_solver;
      }
      else {
         DSmoother coarse_smoother(coarse_A);
         smoothers[level] = &coarse_smoother;
      }

      SparseMatrix P = CreatePixelProlongation(meshes[level], fe_spaces[level], meshes[level+1], fe_spaces[level+1], ess_dofs[level+1]);
      prolongations[level] = &P;
   }

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
