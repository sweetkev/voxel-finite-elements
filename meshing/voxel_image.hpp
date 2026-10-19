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

class VoxelImage 
{
public:

    // Create mesh structure from image file.
    VoxelImage(string file_);
    // Create structure for graph operations from mesh file.
    VoxelImage(const Mesh &mesh_);

    // Returns the mesh read from the image file, or a copy of the input mesh.
    Mesh MakeMesh();

    // Returns map from element numbers to grid cell indices.
    const vector<CellIndex>& GetElementToCell() const { return element_to_cell; }

private:
    string file;
    unique_ptr<Mesh> input_mesh; // Set only when constructed from a mesh

    vector<CellIndex> element_to_cell;

    string GetFileType();
    Mesh ReadPgm();
    Mesh ReadNrrd();
};