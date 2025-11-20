#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <array>
#include <unordered_map>
#include <tuple>

using namespace mfem;

struct Coord
{
   std::array<int,2> coords{};
   Coord() = default;
   Coord(int i, int j) : coords({i, j}) { }
   int &operator[](int i) { return coords[i]; }
   int operator[](int i) const { return coords[i]; }
   bool operator==(const Coord &other) const { return coords[0] == other[0] && coords[1] == other[1]; }
};

template<> struct std::hash<Coord>
{
   std::size_t operator()(const Coord &c) const noexcept
   {
      return c[0]^(c[1] + 0x9e3779b9 + (c[0]<<6) + (c[0]>>2));
   }
};

class PixelMesh
{
public:
   PixelMesh(const PixelImage &image_);

   PixelMesh(std::string pgm_file) : PixelMesh(PixelImage(pgm_file)) { }

   // PixelMesh(std::string pgm_file);
   // PixelMesh(Mesh &mesh_, const std::unordered_map<Coord, int> &coord_to_vertex_, int width_, int height_);

   /** Creates coarse mesh from current mesh */
   PixelMesh CoarsenMesh();

   // /** Creates prolongation map from current mesh to finer mesh */
   // SparseMatrix CreateProlongation(PixelMesh fine_mesh);

   // /** Returns associated mesh */
   Mesh &GetMesh() { return mesh; }

   const Mesh &GetMesh() const { return mesh; }

   // /** Returns associated coordinate pair to vertex map */
   // const std::unordered_map<Coord, int> &GetVertexMap() { return coord_to_vertex; }

   Coord GetElementCoord(int i) const { return element_to_coord.at(i); }

   int GetElementIndex(Coord coord) const
   {
      auto it = coord_to_element.find(coord);
      if (it != coord_to_element.end())
      {
         return it->second;
      }
      else
      {
         return -1;
      }
   }

private:
   PixelImage image;

   Mesh mesh;
   std::unordered_map<Coord, int> coord_to_element;
   std::unordered_map<int, Coord>
   element_to_coord; // TODO: could be std::vetor<Coord>

   int width = -1;
   int height = -1;

   // /** Makes mesh given PixelImage, with elements given by filled-in pixels */
   // std::tuple<std::unordered_map<Coord, int>, Mesh> MakeMesh(const PixelImage &image);

   // /** Returns true when a pixel adjacent to vertex is filled in. Else returns false. */
   // bool AdjacentPixelFilled(int i, int j, const PixelImage &image);

   // /** Adds vertices to mesh at integer indicies, so pixels will be 1x1. Returns map taking coordinates to vertex number. */
   // void AddVertices(const PixelImage &image, Mesh &mesh, int m, int n);

   // /** Adds quads for filled pixels */
   // void AddQuads(const PixelImage &image, Mesh &mesh, int m, int n);

   // /** Returns true if any adjacent coarse pixel should be filled in. Else returns false */
   // bool PixelNearby(int x, int y, const std::unordered_map<Coord, int> &coord_to_fine_vertex);

   // /** Returns unorderd map that takes coordinates to coarse vertex, and adds vertices to coarse mesh */
   // std::unordered_map<Coord, int> AddCoarseVertices(Mesh &coarse_mesh, int m, int n, const std::unordered_map<Coord, int> &coord_to_fine_vertex);

   // /** Adds quads to coarse mesh for filled pixels */
   // void AddCoarseQuads(Mesh &coarse_mesh, int m, int n, const std::unordered_map<Coord, int> &coord_to_fine_vertex);

   // /** Adds the influence of a vertex in the fine mesh to the prolongation matrix P */
   // void AddVertexInfluence(int x, int y, SparseMatrix &P, const std::unordered_map<Coord, int> &fine_coord_to_vertex);
};

SparseMatrix CreatePixelProlongation(const PixelMesh &coarse_mesh,
                                     const FiniteElementSpace &coarse_fes,
                                     const PixelMesh &fine_mesh,
                                     const FiniteElementSpace &fine_fes);
