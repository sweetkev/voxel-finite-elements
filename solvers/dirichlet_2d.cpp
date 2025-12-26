#include "ppm.hpp"
#include "mfem.hpp"
#include "pixel_mesh.hpp"
#include <fstream>
#include <unordered_map>
#include <tuple>

using namespace mfem;
using namespace std;

void AddProlongationBCs(SparseMatrix &P, Array<int> &fine_ess_dofs);

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

   //make mesh from file
   PixelMesh *fine_mesh = new PixelMesh(pgm_file);

   H1_FECollection fec(order, fine_mesh->GetMesh().Dimension());
   FiniteElementSpace *fine_fes = new FiniteElementSpace(&fine_mesh->GetMesh(),
                                                         &fec);

   Array<int> *fine_ess_dofs = new Array<int>;
   fine_fes->GetBoundaryTrueDofs(*fine_ess_dofs);

   BilinearForm fine_a(fine_fes);
   fine_a.AddDomainIntegrator(new DiffusionIntegrator);
   fine_a.Assemble();

   //initial guess
   GridFunction x(fine_fes);
   x = 0.0;

   //right-hand side
   ConstantCoefficient one(1.0);
   LinearForm b(fine_fes);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   SparseMatrix fine_A;
   Vector X, B;
   fine_a.FormLinearSystem(*fine_ess_dofs, x, b, fine_A, X, B);

   // If nlevels was not specified, coarsen until height or width is 2
   if(nlevels < 0) {
      int width = fine_mesh->GetWidth();
      int height = fine_mesh->GetHeight();
      
      // Fnd the number of coarsenings until height or width is 2 (if square, when we have a 2x2 grid)
      nlevels = max(ceil(log2(width)),ceil(log2(height)));
      cout << "new nlevels: " << to_string(nlevels) << "\n";
   }

   //create multigrid operators
   Array<Operator*> operators(nlevels);
   Array<Solver*> smoothers(nlevels);
   Array<Operator*> prolongations(nlevels - 1);
   Array<bool> own_operators(nlevels);
   Array<bool> own_smoothers(nlevels);
   Array<bool> own_prolongations(nlevels - 1);

   own_operators = false;
   own_smoothers = false;
   own_prolongations = false;

   vector<unique_ptr<PixelMesh>> meshes(nlevels);
   vector<unique_ptr<FiniteElementSpace>> fe_spaces(nlevels);
   vector<unique_ptr<BilinearForm>> forms(nlevels);
   vector<unique_ptr<Array<int>>> ess_dofs(nlevels);

   meshes[nlevels-1].reset(fine_mesh);
   fe_spaces[nlevels-1].reset(fine_fes);
   ess_dofs[nlevels-1].reset(fine_ess_dofs);

   operators[nlevels-1] = &fine_A;

   DSmoother fine_smoother(fine_A);
   smoothers[nlevels-1] = &fine_smoother;

   for (int level = nlevels-2; level >= 0; --level)
   {
      PixelMesh *coarse_mesh = new PixelMesh(meshes[level+1]->CoarsenMesh());
      meshes[level].reset(coarse_mesh);

      FiniteElementSpace *coarse_fes = new FiniteElementSpace(
         &coarse_mesh->GetMesh(), &fec);
      fe_spaces[level].reset(coarse_fes);

      Array<int> *coarse_ess_dofs = new Array<int>;
      coarse_fes->GetBoundaryTrueDofs(*coarse_ess_dofs);
      ess_dofs[level].reset(coarse_ess_dofs);

      BilinearForm *coarse_a = new BilinearForm(coarse_fes);
      coarse_a->AddDomainIntegrator(new DiffusionIntegrator);
      coarse_a->Assemble();

      OperatorHandle coarse_A;
      coarse_a->FormSystemMatrix(*coarse_ess_dofs, coarse_A);
      coarse_A.SetOperatorOwner(false);

      operators[level] = coarse_A.Ptr();
      own_operators[level] = true;

      if (level == 0)
      {
         smoothers[0] = new UMFPackSolver(*coarse_A.As<SparseMatrix>());
      }
      else
      {
         smoothers[level] = new DSmoother(*coarse_A.As<SparseMatrix>());
      }
      own_smoothers[level] = true;

      SparseMatrix *P = new SparseMatrix(
         CreatePixelProlongation(*meshes[level], *fe_spaces[level],
                                 *meshes[level+1], *fe_spaces[level+1],
                                 *ess_dofs[level+1]));
      prolongations[level] = P;
      own_prolongations[level] = true;
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
   fine_mesh->GetMesh().Save("mesh.mesh");
}
