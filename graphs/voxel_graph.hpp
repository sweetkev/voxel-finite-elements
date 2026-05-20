#pragma once

#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>

using namespace mfem;
using namespace std;

class VoxelGraph
{
public:
    // Create a graph where each vertex of the graph corresponds to an occupied
    // pixel/voxel in the image, and edges are created using vertex adjacency
    // of the pixels/voxels. 
    VoxelGraph(const FiniteElementSpace &fes, const PixelImage &image_);

    // Returns the graph over a coarsened mesh. Each grid cell of the coarse
    // mesh is obtained by coarsening a 2x2(x2) block of grid cells in the fine
    // mesh. Each connected component of the graph over a coarse grid cell 
    // becomes a new vertex in the coarse graph.
    VoxelGraph CoarsenGraph();

    const vector<set<int>> &GetDofGroups() const { return dof_groups; }

    const vector<int> &GetBrokenToTrueDofMap() const 
    { 
        return broken_to_true_dof; 
    }

    const vector<int> &GetFineToCoarseElementMap() const 
    { 
        return fine_to_coarse_element_map;
    }

    const Coord &GetElementCoord(int e) const { return element_to_coord[e]; }

    int Size() const { return graph.Size(); }

private:
    // Graph structure
    Table graph;
    vector<Coord> grid_cells; // Occupied grid cells

    // Maps between coordinates and elements
    unordered_map<Coord, vector<int>> coord_to_element;
    vector<Coord> element_to_coord;

    // Fine-to-coarse element map. This is the "labeling" of each fine element
    // that determines which coarse element it belongs to.
    vector<int> fine_to_coarse_element_map;

    // Reference FES needed to determine DoF identification
    FiniteElementSpace reference_fes;

    // Vector with size equal to number of true dofs whose entries contain the
    // sets of broken dofs identified with the true dof.    
    vector<set<int>> dof_groups;

    // Maps each broken dof index to the true dof index.
    // Size is (ne * dofs_per_element).
    std::vector<int> broken_to_true_dof;

    // Matrix whose (i,j) entry denotes the local dof number on reference
    // element j which is identified with the local dof i on
    // the central reference element.
    DenseMatrix local_to_neighbor_dof_map;    

    VoxelGraph(Table &graph, vector<Coord> &grid_cells,
               unordered_map<Coord, vector<int>> &coord_to_element,
               vector<Coord> &element_to_coord, 
               const FiniteElementSpace &reference_fes);

    // Creates the matrix whose ith row represents a local dof number
    // and jth column represents a reference element in the reference mesh.
    void CreateReferenceDofMap();
    
    // Creates the vector with size equal to number of true dofs whose entries
    // contain the sets of broken dofs identified with the true dof.
    void CreateDofGroups();

    // Returns the local dof number on neighbor element that corresponds
    // to the local dof number e_dof on element e. If the dof is not shared
    // between elements, returns -1
    int GetNeighborDof(int e, int e_dof, int neighbor);
};