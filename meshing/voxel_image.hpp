#pragma once

#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <array>
#include <memory>
#include <unordered_map>
#include <tuple>

using namespace mfem;
using namespace std;

struct CellIndex
{
   std::array<int,3> ijk{};   // k == 0 for 2D cells
   CellIndex() = default;
   CellIndex(int i, int j, int k = 0) : ijk({i, j, k}) { }
   int &operator[](int i) { return ijk[i]; }
   int operator[](int i) const { return ijk[i]; }
   bool operator==(const CellIndex &other) const { return ijk == other.ijk; }
   bool operator!=(const CellIndex &other) const { return ijk != other.ijk; }
};

template<> struct std::hash<CellIndex>
{
    std::size_t operator()(const CellIndex &c) const noexcept
    {
        std::size_t h = 0;
        for (int x : c.ijk) { 
            h ^= std::hash<int>{}(x) + 0x9e3779b9 + (h << 6) + (h >> 2); 
        }
        return h;
    }
};

// Returns the grid cell containing each element of the mesh, computed from the
// element's vertex coordinates. The grid cell (i,j,k) is the element whose 
// lower left corner is nearest to h*(i,j,k).
vector<CellIndex> ComputeElementCells(const Mesh &mesh, real_t h);

class VoxelImage 
{
public:

    // Create mesh structure from image file.
    VoxelImage(string file_);
    // Create structure for graph operations from mesh file.
    VoxelImage(const Mesh &mesh_);

    // Returns the mesh read from the image file, or a copy of the input mesh.
    Mesh MakeMesh();

    // Returns the side length of elements.
    real_t GetCellSize() const
    {
        return h;
    }

private:
    string file;
    unique_ptr<Mesh> input_mesh; // Set only when constructed from a mesh

    real_t h = 0.0; // Set by the constructor or MakeMesh()

    string GetFileType();
    Mesh ReadPgm();
    Mesh ReadNrrd();
};