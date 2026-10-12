#pragma once

#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <array>
#include <unordered_map>
#include <tuple>

using namespace mfem;

struct CellIndex
{
   std::array<int,2> ijk{};
   CellIndex() = default;
   CellIndex(int i, int j) : ijk({i, j}) { }
   int &operator[](int i) { return ijk[i]; }
   int operator[](int i) const { return ijk[i]; }
   bool operator==(const CellIndex &other) const
      { return ijk[0] == other[0] && ijk[1] == other[1]; }
   bool operator!=(const CellIndex &other) const
      { return ijk[0] != other[0] || ijk[1] != other[1]; }
};

template<> struct std::hash<CellIndex>
{
   std::size_t operator()(const CellIndex &c) const noexcept
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

   /** Returns dimension of mesh */
   int Dimension() const { return mesh.Dimension(); }

   /** Returns mesh width */
   const int GetWidth() const { return width; }

   /** Returns mesh height */
   const int GetHeight() const { return height; }

   /** Returns the grid cell of element i */
   CellIndex GetElementCell(int i) const { return element_to_cell[i]; }

   /**
    * Returns the index of the element in the given grid cell
    */
   int GetElementIndex(CellIndex cell) const
   {
      auto it = cell_to_element.find(cell);
      if (it != cell_to_element.end())
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
   std::unordered_map<CellIndex, int> cell_to_element;
   std::vector<CellIndex> element_to_cell;

   int width = -1;
   int height = -1;
};

SparseMatrix CreatePixelProlongation(const PixelMesh &coarse_mesh,
                                     const FiniteElementSpace &coarse_fes,
                                     const PixelMesh &fine_mesh,
                                     const FiniteElementSpace &fine_fes,
                                     const Array<int> &fine_ess_dofs);

void EnforceProlongationBCs(SparseMatrix &P, const Array<int> &fine_ess_dofs);
