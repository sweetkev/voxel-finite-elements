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

   /** Returns coarsened mesh obtained from current mesh */
   PixelMesh CoarsenMesh() { return PixelMesh(image.Coarsen()); }


   /** Returns associated mesh */
   Mesh &GetMesh() { return mesh; }

   /** Returns associated mesh as const */
   const Mesh &GetMesh() const { return mesh; }

   /** Returns mesh width */
   const int GetWidth() const { return width; }

   /** Returns mesh height */
   const int GetHeight() const { return height; }

   /** Returns the bottom-left coordinates of element i */
   Coord GetElementCoord(int i) const { return element_to_coord[i]; }

   /** 
    * Returns the element index of the element with bottom-left 
    *  coordinates coord 
    */
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
   std::vector<Coord> element_to_coord;

   int width = -1;
   int height = -1;
};

SparseMatrix CreatePixelProlongation(const PixelMesh &coarse_mesh,
                                     const FiniteElementSpace &coarse_fes,
                                     const PixelMesh &fine_mesh,
                                     const FiniteElementSpace &fine_fes,
                                     const Array<int> &fine_ess_dofs);

void EnforceProlongationBCs(SparseMatrix &P, const Array<int> &fine_ess_dofs);
