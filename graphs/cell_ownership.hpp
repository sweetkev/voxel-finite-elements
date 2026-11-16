#pragma once

#include "mfem.hpp"
#include "voxel_image.hpp"
#include <unordered_map>
#include <vector>

using namespace mfem;
using namespace std;

// Returns the grid cell on coarse_level containing grid cell `cell` on
// current_level.
CellIndex AncestorCell(const CellIndex &cell, int current_level,
                       int coarse_level);

// Determines, for every coarse level, which MPI rank owns each coarse grid
// cell. The owner of a coarse cell is the lowest rank holding any finest-level
// element inside it; it is the rank that computes the connected components
// (coarse elements) on that cell.
//
//   1. Each rank sends (level, coarse cell) for every distinct ancestor of its
//      fine cells to a directory rank. Currently the root is the directory
//      rank for every cell.
//   2. The directory rank determines the owner of each coarse cell, which is
//      the smallest-indexed rank with fine elements on that cell.
//   3. The directory replies to every contributor with the owner, and to the
//      owner additionally with the full list of contributors.
class CellOwnership
{
public:
    CellOwnership(MPI_Comm comm, const vector<CellIndex> &fine_cells,
                  int nlevels);

    int NumLevels() const { return nlevels; }

    // Owner of grid cell `cell` on coarse `level`. Only valid for cells
    // containing local fine elements.
    int GetOwner(int level, const CellIndex &cell) const
    {
        return owner[level].at(cell);
    }

    // Cells on coarse `level` that contain local fine elements, mapped to their
    // owner.
    const unordered_map<CellIndex, int> &GetOwners(int level) const
    {
        return owner[level];
    }

    // Cells on coarse `level` owned by this rank, mapped to the sorted list of
    // ranks holding fine elements inside them (including this rank).
    const unordered_map<CellIndex, vector<int>> &GetOwnedCells(int level) const
    {
        return contributors[level];
    }

private:
    MPI_Comm comm;
    int nlevels;
    // Indexed by coarse level 0, ..., nlevels - 2
    vector<unordered_map<CellIndex, int>> owner;
    vector<unordered_map<CellIndex, vector<int>>> contributors;
};
