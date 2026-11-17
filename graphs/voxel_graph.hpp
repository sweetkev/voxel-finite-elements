#pragma once

#include "voxel_image.hpp"
#include "cell_ownership.hpp"
#include "mfem.hpp"
#include <unordered_map>

using namespace mfem;
using namespace std;

class VoxelGraph
{
public:
    // Create a graph where each vertex of the graph corresponds to an occupied
    // pixel/voxel in the image, and edges are created using vertex adjacency
    // of the pixels/voxels. element_to_cell gives the grid cell of each
    // element of pfes's mesh. Each rank holds the elements of its part of the 
    // mesh, and neighbors on other ranks are found through the shared mesh 
    // vertices.
    VoxelGraph(const ParFiniteElementSpace &pfes,
               vector<CellIndex> element_to_cell,
               const FiniteElementSpace &reference_fes);

    // Returns the graph over a coarsened mesh. Each grid cell of the coarse 
    // mesh is obtained by coarsening a 2x2(x2) block of grid cells in the finer
    // mesh. Each connected component of the graph over a coarse grid cell 
    // becomes a new vertex in the coarse graph.
    //
    // All of the coarsening work for a coarse grid cell (connected components
    // and coarse adjacency) is done on the cell's owner from ownership, and
    // the resulting coarse elements live on that rank:
    //   1. Every fine element is sent with its neighbor list to the owner of
    //      its coarse grid cell.
    //   2. Each owner labels the connected components of its cells, and tells
    //      the rank holding each fine element which coarse element it joined.
    //   3. Each owner asks the ranks holding neighbors outside its cells for
    //      their coarse elements, then builds the coarse adjacency.
    VoxelGraph CoarsenGraph(const CellOwnership &ownership, int coarse_level);

    // Global number of local element e. Elements are numbered contiguously by
    // rank (like the values in CSR format).
    int GlobalElementNumber(int e) const { return element_offsets[rank] + e; }

    // Rank holding the element with global number global_e.
    int ElementRank(int global_e) const;

    int GlobalNE() const { return element_offsets.back(); }

    // Adjacency between local elements, in local element numbering.
    const Table &GetGraph() const { return graph; }

    // Global numbers of the neighbors of each local element that live on
    // other ranks (one row per local element).
    const Table &GetOffRankNeighbors() const
    {
        return off_rank_neighbors;
    }

    // Number of local dofs: the dofs on this rank after identifying the broken
    // dofs shared by neighboring local elements. A dof shared with an element
    // on another rank has a local copy on each rank, and GetDofTrueDofMatrix
    // maps the copies to one global true dof.
    int NumLocalDofs() const { return num_local_dofs; }

    // Maps each broken dof to its local dof.
    const vector<int> &GetBrokenToLocalDofMap() const
    {
        return broken_to_local_dof;
    }

    // Boolean matrix mapping global true dofs to local dofs. On the finest 
    // level it is the matrix of the ParFiniteElementSpace passed to the 
    // constructor.
    HypreParMatrix &GetDofTrueDofMatrix() const { return *dof_truedof; }

    // Global number of the coarse element containing each local element.
    const vector<int> &GetFineToCoarseElementMap() const 
    { 
        return fine_to_coarse_element_map;
    }

    // Global numbers of the fine elements contained in each coarse element
    // owned by this rank (local coarse numbering).
    const vector<set<int>> &GetCoarseToFineElementsMap() const
    {
        return coarse_to_fine_elements_map;
    }

    const CellIndex &GetElementCell(int e) const { return element_to_cell[e]; }

    int Size() const { return graph.Size(); }

private:
    MPI_Comm comm;
    int rank;

    // element_offsets[r] is the global number of the first element on rank r;
    // size is nranks + 1.
    vector<int> element_offsets;

    // Graph structure. graph holds edges between local elements;
    // off_rank_neighbors holds the global numbers of neighbors on other ranks.
    Table graph;
    Table off_rank_neighbors;
    vector<CellIndex> grid_cells; // Occupied grid cells

    // Maps between grid cells and elements
    unordered_map<CellIndex, vector<int>> cell_to_elements;
    vector<CellIndex> element_to_cell;

    // Fine-to-coarse element map. This is the "labeling" of each fine element
    // that determines which coarse element it belongs to (global coarse
    // numbering).
    vector<int> fine_to_coarse_element_map;
    vector<set<int>> coarse_to_fine_elements_map;

    // Reference FES needed to determine DoF identification
    const FiniteElementSpace &reference_fes;

    // Number of local dofs (see NumLocalDofs).
    int num_local_dofs = 0;

    // Maps each broken dof index to the local dof index.
    // Size is (ne * dofs_per_element).
    vector<int> broken_to_local_dof;

    // Local dof to global true dof matrix. Points to owned_dof_truedof on
    // coarse levels and to the ParFiniteElementSpace's matrix on the finest
    // level.
    HypreParMatrix *dof_truedof = nullptr;
    unique_ptr<HypreParMatrix> owned_dof_truedof;

    // Matrix whose (i,j) entry denotes the local dof number on reference
    // element j which is identified with the local dof i on
    // the central reference element.
    DenseMatrix local_to_neighbor_dof_map;    

    VoxelGraph(MPI_Comm comm, vector<int> element_offsets, Table graph,
               Table off_rank_neighbors,
               vector<CellIndex> grid_cells,
               unordered_map<CellIndex, vector<int>> cell_to_elements,
               vector<CellIndex> element_to_cell, 
               const FiniteElementSpace &reference_fes);

    // Finds off_rank_neighbors: elements on other ranks sharing a mesh vertex
    // with a local element. Collective.
    void FindOffRankNeighbors(const ParFiniteElementSpace &pfes);

    // Creates the matrix whose ith row represents a local dof number
    // and jth column represents a reference element in the reference mesh.
    void CreateReferenceDofMap();
    
    // Identifies broken dofs shared between neighboring local elements, then
    // builds broken_to_local_dof and sets num_local_dofs.
    void CreateBrokenToLocalDofMap();

    // Builds owned_dof_truedof on a coarse level by identifying local dofs
    // with the matching dofs of neighboring elements on other ranks.
    //   1. Each local dof is labeled with a global number given by its local 
    //      number plus the number of local dofs on lower ranks.
    //   2. Each rank sends the labels of local dofs on every element with 
    //      off-rank neighbors to the neighbors' ranks. Each received label 
    //      replaces the label of the matching local dof if it is smaller. This 
    //      repeats until no label changes, so every copy of a dof ends with the
    //      smallest label among its copies.
    //   3. The rank holding the local dof with that smallest label owns the
    //      true dof. Owned true dofs are numbered globally, and the other ranks
    //      ask the owners for the true dof numbers of their copies.
    void CreateDofTrueDofMatrix();

    // Returns the local dof number on neighbor element that corresponds
    // to the local dof number e_dof on element e. If the dof is not shared
    // between elements, returns -1
    int GetNeighborDof(int e, int e_dof, int neighbor) const
    {
        return GetNeighborDof(element_to_cell[e], e_dof,
                              element_to_cell[neighbor]);
    }

    // Same as above, for elements in grid cells e_cell and neighbor_cell.
    int GetNeighborDof(const CellIndex &e_cell, int e_dof,
                       const CellIndex &neighbor_cell) const;

    // Reindex DOFs to match MFEM's labeling.
    void ReindexDofs(const FiniteElementSpace &fes);
};