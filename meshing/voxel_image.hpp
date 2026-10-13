#pragma once

#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <array>
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

    VoxelImage(string file_);
    // TODO: VoxelImage (Mesh mesh);

    Mesh MakeMesh();

    const vector<int>& GetCellDims() const { return cell_dims; }

    // Returns true if cell (i,j,k) is shaded. k=0 for 2D images
    bool operator()(int i, int j, int k=0) const;
    // Returns true if cell i is shaded. 
    bool operator[](int i) const;

private:
    string file;

    vector<int> cell_dims;

    vector<bool> data;
    
    string GetFileType();
    Mesh ReadPgm();
    Mesh ReadNrrd();
};