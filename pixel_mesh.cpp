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
         }
      }
   }

   the_mesh.RemoveUnusedVertices();
   the_mesh.FinalizeMesh(); // The boundary elements will be generated here.

   mesh.Swap(the_mesh, true);
}

PixelMesh PixelMesh::CoarsenMesh()
{
   return PixelMesh(image.Coarsen());

   //    int m = width, n = height;
   //    std::unordered_map<Coord, int> coord_to_fine_vertex = coord_to_vertex;

   //    //dimension of domain and ambient space
   //    int dim = 2, sdim = 2;

   //    //number of vertices, elements, and boundary elements. Will be allocated by FinalizeMesh() later
   //    int nv = 0;
   //    int ne = 0;
   //    int nb = 0;

   //    //initialize mesh
   //    Mesh coarse_mesh(dim, nv, ne, nb, sdim);

   //    //add vertices
   //    std::unordered_map<Coord, int> coord_to_coarse_vertex = AddCoarseVertices(
   //                                                               coarse_mesh, m, n, coord_to_fine_vertex);

   //    //add quads
   //    AddCoarseQuads(coarse_mesh, m, n, coord_to_fine_vertex, coord_to_coarse_vertex);

   //    coarse_mesh.FinalizeMesh();
   //    coarse_mesh.Save("coarse_mesh.mesh");

   //    PixelMesh c_mesh(coarse_mesh, coord_to_coarse_vertex, m, n);
   //    return c_mesh;
}

// PixelMesh::PixelMesh(Mesh &mesh_, const std::unordered_map<Coord, int> &coord_to_vertex_, int width_, int height_)
//     : mesh(mesh_), coord_to_vertex(coord_to_vertex_), width(width_), height(height_)
// {}



// SparseMatrix PixelMesh::CreateProlongation(PixelMesh fine_mesh) {
//     std::unordered_map<Coord, int> fine_coord_to_vertex = fine_mesh.GetVertexMap();

//     //initialize matrix
//     int nrows = fine_coord_to_vertex.size();
//     int ncols = coord_to_vertex.size();
//     SparseMatrix P(nrows,ncols);

//     //add entries
//     int p = std::ceil(0.5*(width))+1, q = std::ceil(0.5*(height))+1;

//     for(int j = 0; j < q; j++) {
//         for(int i = 0; i < p; i++) {
//             int x = 2*i, y = 2*(q-j-1);
//             std::string coord = Coord(x, y);
//             if(coord_to_vertex.count(coord) != 0) {
//                 AddVertexInfluence(x,y,P,fine_coord_to_vertex);
//             }
//         }
//     }

//     //finalize and return
//     P.Finalize();
//     return P;
// }

// std::tuple<std::unordered_map<Coord, int>, Mesh> PixelMesh::MakeMesh(PixelImage image) {
//     //dimension of domain and ambient space
//     int dim = 2, sdim = 2;

//     //number of vertices, elements, and boundary elements. Will be allocated by FinalizeMesh() later
//     int nv = 0;
//     int ne = 0;
//     int nb = 0;

//     //initialize mesh
//     Mesh fine_mesh(dim, nv, ne, nb, sdim);

//     //add vertices
//     int m = image.Width(), n = image.Height();
//     std::unordered_map<Coord, int> coord_to_vertex = AddVertices(image, fine_mesh, m, n);

//     //add quads
//     AddQuads(image, fine_mesh, m, n, coord_to_vertex);

//     fine_mesh.FinalizeMesh();
//     fine_mesh.Save("fine_mesh.mesh");

//     return std::make_tuple(coord_to_vertex, fine_mesh);
// }

// bool PixelMesh::AdjacentPixelFilled(int i, int j, PixelImage image) {
//     int m = image.Width(), n = image.Height();

//     //bottom row cases
//     if(j == n) {
//         //bottom left vertex
//         if(i == 0) {
//             return image.operator()(i,j-1) != 0;
//         }
//         //bottom right vertex
//         else if(i == m) {
//             return image.operator()(i-1,j-1) != 0;
//         }
//         //other bottom vertices
//         else {
//             return (image.operator()(i-1,j-1) != 0) || (image.operator()(i,j-1) != 0);
//         }
//     }

//     //top row cases
//     if(j == 0) {
//         //top left vertex
//         if(i == 0) {
//             return image.operator()(i,j) != 0;
//         }
//         //top right vertex
//         else if(i == m) {
//             return image.operator()(i-1,j) != 0;
//         }
//         //other top vertices
//         else {
//             return (image.operator()(i,j) != 0) || (image.operator()(i-1,j) != 0);
//         }
//     }

//     //all other rows
//     //left vertices
//     if(i == 0) {
//         return (image.operator()(i,j-1) != 0) || (image.operator()(i,j) != 0);
//     }
//     //right vertices
//     else if(i == m) {
//         return (image.operator()(i-1,j-1) != 0) || (image.operator()(i-1,j) != 0);
//     }
//     //all other vertices (interior vertices)
//     else {
//         return (image.operator()(i-1,j-1) != 0) || (image.operator()(i-1,j) != 0) || (image.operator()(i,j-1) != 0) || (image.operator()(i,j) != 0);
//     }
// }

// void PixelMesh::AddVertices(PixelImage image, Mesh &mesh, int m, int n) {
//     int vertex = 0;
//     for(int j = 0; j < n+1; j++) {
//         for(int i = 0; i < m+1; i++) {
//             if(AdjacentPixelFilled(i,j,image)) {
//                 int x = i, y = n-j;
//                 mesh.AddVertex(x,y);
//                 std::string coord = Coord(x, y);
//                 coord_to_vertex.insert({coord, vertex});
//                 vertex++;
//             }
//         }
//     }
// }

// void PixelMesh::AddQuads(PixelImage image, Mesh &mesh, int m, int n) {
//     for(int j = 0; j < n; j++) {
//         for(int i = 0; i < m; i++) {
//             //add quad for filled in pixels
//             if(image.operator()(i,j) != 0) {
//                 int x = i, y = n-j;
//                 int v1 = coord_to_vertex[Coord(x, y)];
//                 int v2 = coord_to_vertex[Coord(x, y-1)];
//                 int v3 = coord_to_vertex[Coord(x+1, y-1)];
//                 int v4 = coord_to_vertex[Coord(x+1, y)];
//                 mesh.AddQuad(v1,v2,v3,v4);
//             }
//         }
//     }
// }

// bool PixelMesh::PixelNearby(int x, int y, std::unordered_map<Coord, int> coord_to_fine_vertex) {
//     Coord coord1(x - 1, y - 1);
//     Coord coord2(x + 1, y - 1);
//     Coord coord3(x + 1, y + 1);
//     Coord coord4(x - 1, y + 1);
//     if((coord_to_fine_vertex.count(coord1) != 0) || (coord_to_fine_vertex.count(coord2) != 0) || (coord_to_fine_vertex.count(coord3) != 0) || (coord_to_fine_vertex.count(coord4) != 0)) {
//         return true;
//     }

//     return false;
// }

// std::unordered_map<Coord, int> PixelMesh::AddCoarseVertices(Mesh &coarse_mesh, int m, int n, std::unordered_map<Coord, int> coord_to_fine_vertex) {
//     int p = std::ceil(0.5*(m))+1, q = std::ceil(0.5*(n))+1;

//     std::unordered_map<Coord, int> coord_to_coarse_vertex;

//     int vertex = 0;
//     for(int j = 0; j < q; j++) {
//         for(int i = 0; i < p; i++) {
//             int x = 2*i, y = 2*(q-j-1);
//             std::string coord = Coord(x, y);
//             if(PixelNearby(x,y,coord_to_fine_vertex)) {
//                 coarse_mesh.AddVertex(x,y);
//                 coord_to_coarse_vertex.insert({coord,vertex});
//                 vertex++;
//             }
//         }
//     }

//     return coord_to_coarse_vertex;
// }

// void PixelMesh::AddCoarseQuads(Mesh &coarse_mesh, int m, int n, const std::unordered_map<Coord, int> &coord_to_fine_vertex) {
//     int p = std::ceil(0.5*m), q = std::ceil(0.5*n);

//     for(int j = 0; j < q; j++) {
//         for(int i = 0; i < p; i++) {
//             //add quad for filled in pixels
//             int x = 2*i, y = 2*(q-j-1);
//             Coord coord(x+1, y-1);
//             if(coord_to_fine_vertex.count(coord) != 0) {
//                 int v1 = coord_to_vertex[Coord(x, y)];
//                 int v2 = coord_to_vertex[Coord(x, y-2)];
//                 int v3 = coord_to_vertex[Coord(x+2, y-2)];
//                 int v4 = coord_to_vertex[Coord(x+2, y)];
//                 coarse_mesh.AddQuad(v1,v2,v3,v4);
//             }
//         }
//     }
// }

// void PixelMesh::AddVertexInfluence(int x, int y, SparseMatrix &P, const std::unordered_map<Coord, int> &fine_coord_to_vertex) {
//     Coord coarse_coord(x, y);
//     int i,j;
//     j = coord_to_vertex[coarse_coord];

//     Coord coord;

//     //current (middle middle) vertex
//     coord = coarse_coord;
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,1);
//     }

//     //bottom left vertex
//     coord = Coord(x-1, y-1);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.25);
//     }

//     //bottom middle vertex
//     coord = Coord(x, y-1);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.5);
//     }

//     //bottom right vertex
//     coord = Coord(x+1, y-1);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.25);
//     }

//     //middle left vertex
//     coord = Coord(x-1, y);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.5);
//     }

//     //middle right vertex
//     coord = Coord(x+1, y);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.5);
//     }

//     //top left vertex
//     coord = Coord(x-1, y+1);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.25);
//     }

//     //top middle vertex
//     coord = Coord(x, y+1);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.5);
//     }

//     //top right vertex
//     coord = Coord(x+1, y+1);
//     if(fine_coord_to_vertex.count(coord) != 0) {
//         i = fine_coord_to_vertex[coord];
//         P.Add(i,j,0.25);
//     }
// }

struct ChildIndex
{
   int fine_element;
   int sub_element_position;
};

SparseMatrix CreatePixelProlongation(const PixelMesh &coarse_mesh,
                                     const FiniteElementSpace &coarse_fes,
                                     const PixelMesh &fine_mesh,
                                     const FiniteElementSpace &fine_fes)
{
   SparseMatrix P(fine_fes.GetTrueVSize(), coarse_fes.GetTrueVSize());

   const int coarse_ne = coarse_mesh.GetMesh().GetNE();

   std::vector<ChildIndex> children(fine_mesh.GetMesh().GetNE());
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
   //    static const double A = 0.0, B = 0.5, C = 1.0;
   //    // NOTE: as opposed to Mesh::UniformRefinement2D_base, these are ordered
   //    // lexicographically
   //    static double quad_children[2*4*4] =
   //    {
   //       A,A, B,A, B,B, A,B, // lower-left
   //       B,A, C,A, C,B, B,B, // lower-right
   //       A,B, B,B, B,C, A,C,  // upper-left
   //       B,B, C,B, C,C, B,C // upper-right
   //    };

   //    const FiniteElement &fe = *coarse_fes.GetFE(0);
   //    const int n_loc_dof = fe.GetDof();

   //    IsoparametricTransformation isotr;
   //    isotr.SetIdentityTransformation(fe.GetGeomType());

   //    local_P.SetSize(n_loc_dof, n_loc_dof, 4);
   //    local_R.SetSize(n_loc_dof, n_loc_dof, 4);
   //    for (int i = 0; i < 4; ++i)
   //    {
   //       DenseMatrix pmat(quad_children + i*2*4, 2, 4);
   //       isotr.SetPointMat(pmat);
   //       fe.GetLocalInterpolation(isotr, local_P(i));
   //       fe.GetLocalRestriction(isotr, local_R(i));
   //    }

}
