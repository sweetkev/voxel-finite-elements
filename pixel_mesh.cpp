#include "ppm.hpp"
#include "mfem.hpp"
#include "pixel_mesh.hpp"
#include <iostream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

PixelMesh::PixelMesh(const PixelImage &image_) : image(image_)
{
   width = image.Width();
   height = image.Height();

   //dimension of domain and ambient space
   int dim = 2, sdim = 2;

   //number of vertices, elements, and boundary elements. Will be allocated by FinalizeMesh() later
   int nv = (width + 1)*(height + 1);
   int ne = 0;
   for (int i = 0; i < width*height; ++i)
   {
      if (image[i] != 0) { ++ne; }
   }
   int nb = 0;

   //initialize mesh
   Mesh the_mesh(dim, nv, ne, nb, sdim);

   // Add all verties in the grid, even if they are unused. We will remove
   // unused vertices later.
   real_t h = 1.0 / width;
   for (int j = 0; j < height + 1; ++j)
   {
      for (int i = 0; i < width + 1; ++i)
      {
         the_mesh.AddVertex(i*h, j*h);
      }
   }

   // Add the elements.
   element_to_coord.resize(width*height);
   int e = 0;
   for (int j = 0; j < height; ++j)
   {
      for (int i = 0; i < width; ++i)
      {
         if (image(i, j) != 0)
         {
            Coord coord(i, j);
            coord_to_element[coord] = e;
            element_to_coord[e] = coord;

            const int v1 = j*(width + 1) + i;
            const int v2 = v1 + 1;
            const int v3 = v1 + width + 1 + 1;
            const int v4 = v3 - 1;
            the_mesh.AddQuad(v1, v2, v3, v4);
            e++;
         }
      }
   }

   the_mesh.RemoveUnusedVertices();
   the_mesh.FinalizeMesh(); // The boundary elements will be generated here.

   mesh.Swap(the_mesh, true);
}

/* Each child index corresponds to a fine element. The sub_element_position indicates
   the position of the fine element in realtion to its parental coarse element.
   sub_element_position is lexicographic left-to-right, bottom-to-top, i.e.
   2 3
   0 1
*/
struct ChildIndex
{
   int fine_element;
   int sub_element_position;
};

/** Returns Prolongation Matrix from a coarse mesh to a finer mesh */
SparseMatrix CreatePixelProlongation(const PixelMesh &coarse_mesh,
                                     const FiniteElementSpace &coarse_fes,
                                     const PixelMesh &fine_mesh,
                                     const FiniteElementSpace &fine_fes,
                                     const Array<int> &fine_ess_dofs)
{
   // Initialize prolongation matrix

   int nrows = fine_fes.GetTrueVSize(), ncols = coarse_fes.GetTrueVSize();
   SparseMatrix P(nrows,ncols);

   // Organize structure between coarse element "parents" and fine element "children"

   const int coarse_ne = coarse_mesh.GetMesh().GetNE();
   const int fine_ne = fine_mesh.GetMesh().GetNE();

   /* children and children_offsets utilize a CSR-like structure. For a
      coarse element i, its children can be accessed by the indicies from
      children_offsets[i] up to (not including) children_offsets[i+1].
   */


   std::vector<ChildIndex> children(fine_ne);
   std::vector<int> children_offsets(coarse_ne + 1);

   int offset = 0;
   for (int i = 0; i < coarse_ne; ++i)
   {
      children_offsets[i] = offset;

      const Coord coarse_coord = coarse_mesh.GetElementCoord(i);
      Coord fine_coord(2*coarse_coord[0], 2*coarse_coord[1]);

      for (int jj = 0; jj < 2; ++jj)
      {
         fine_coord[1] += jj;
         for (int ii = 0; ii < 2; ++ii)
         {
            fine_coord[0] += ii;
            const int fine_idx = fine_mesh.GetElementIndex(fine_coord);
            if (fine_idx >= 0)
            {
               children[offset] = {fine_idx, ii + 2*jj};
               ++offset;
            }
            fine_coord[0] -= ii;
         }
         fine_coord[1] -= jj;
      }
   }
   children_offsets.back() = offset;

   // Set up the local prolongation matrices

   // See Mesh::UniformRefinement2D_base
   static const double A = 0.0, B = 0.5, C = 1.0;
   // NOTE: as opposed to Mesh::UniformRefinement2D_base, these are ordered
   // lexicographically
   static double quad_children[2*4*4] =
   {
      A,A, B,A, B,B, A,B, // lower-left
      B,A, C,A, C,B, B,B, // lower-right
      A,B, B,B, B,C, A,C,  // upper-left
      B,B, C,B, C,C, B,C // upper-right
   };

   const FiniteElement &fe = *coarse_fes.GetFE(0);
   const int n_loc_dof = fe.GetDof();

   IsoparametricTransformation isotr;
   isotr.SetIdentityTransformation(fe.GetGeomType());

   DenseTensor local_P; // local prolongation
   // DenseTensor local_R; // local restriction
   local_P.SetSize(n_loc_dof, n_loc_dof, 4);
   // local_R.SetSize(n_loc_dof, n_loc_dof, 4);
   for (int i = 0; i < 4; ++i)
   {
      DenseMatrix pmat(quad_children + i*2*4, 2, 4);
      isotr.SetPointMat(pmat);
      fe.GetLocalInterpolation(isotr, local_P(i));
      // fe.GetLocalRestriction(isotr, local_R(i));
   }

   for (int i = 0; i < coarse_ne; ++i)
   {
      for (int j = children_offsets[i]; j < children_offsets[i+1]; ++j)
      {
         ChildIndex child = children[j];
         int child_position = child.sub_element_position;
         int child_element = child.fine_element;

         Array<int> rows, cols;
         fine_fes.GetElementDofs(child_element, rows);
         coarse_fes.GetElementDofs(i, cols);

         P.SetSubMatrix(rows,cols,local_P(child_position));
      }
   }

   EnforceProlongationBCs(P,fine_ess_dofs);

   return P;

}

/** Removes rows from prolongation P that correspond with boundary dofs */
void EnforceProlongationBCs(SparseMatrix &P, const Array<int> &fine_ess_dofs)
{
   int nbdofs = fine_ess_dofs.Size();

   for (int i=0; i < nbdofs; i++)
   {
      P.EliminateRow(fine_ess_dofs[i]);
   }
}
