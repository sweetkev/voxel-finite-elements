#pragma once

#include "mfem.hpp"
#include "voxel_image.hpp"
#include <unordered_map>
#include <vector>

using namespace mfem;
using namespace std;

// All-to-all exchange: send[r] goes to rank r. Returns the
// concatenation of all received buffers ordered by source rank, with the number
// of values received from each rank in recv_counts.
template <typename T>
vector<T> Exchange(MPI_Comm comm, const vector<vector<T>> &send,
                   vector<int> &recv_counts)
{
    const MPI_Datatype type = MPITypeMap<T>::mpi_type;
    const int nranks = send.size();
    vector<int> send_counts(nranks), send_displs(nranks + 1, 0);
    for (int r = 0; r < nranks; ++r)
    {
        send_counts[r] = send[r].size();
        send_displs[r + 1] = send_displs[r] + send_counts[r];
    }

    recv_counts.assign(nranks, 0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                 recv_counts.data(), 1, MPI_INT, comm);

    vector<int> recv_displs(nranks + 1, 0);
    for (int r = 0; r < nranks; ++r)
    {
        recv_displs[r + 1] = recv_displs[r] + recv_counts[r];
    }

    vector<T> send_buf(send_displs[nranks]);
    for (int r = 0; r < nranks; ++r)
    {
        copy(send[r].begin(), send[r].end(), send_buf.begin() + send_displs[r]);
    }

    vector<T> recv_buf(recv_displs[nranks]);
    MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(), type,
                  recv_buf.data(), recv_counts.data(), recv_displs.data(), type,
                  comm);
    return recv_buf;
}

// Returns the grid cell on coarse_level containing grid cell `cell` on
// current_level.
CellIndex AncestorCell(const CellIndex &cell, int current_level,
                       int coarse_level);

// Position of `cell` among the children of its parent cell one level coarser,
// in the ordering of MFEM's uniform refinement of a quadrilateral/hexahedron
// (the order of CoarseFineTransformations::point_matrices).
int ChildIndex(const CellIndex &cell);

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
