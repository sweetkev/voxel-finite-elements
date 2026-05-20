# pragma once

#include "mfem.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"

using namespace mfem;
using namespace std;

class VoxelGraphHierarchy
{
public:
    VoxelGraphHierarchy(unique_ptr<VoxelGraph> fine_graph, int nlevels,
                        const FiniteElementSpace &reference_fes, real_t h);

    vector<unique_ptr<VoxelGraphOperator>> &GetGraphOperators() { return graph_operators; }
    vector<unique_ptr<SparseMatrix>> &GetProlongations() { return prolongations; }

private:
    vector<unique_ptr<VoxelGraph>> graphs;
    vector<unique_ptr<VoxelGraphOperator>> graph_operators;
    vector<unique_ptr<SparseMatrix>> prolongations;

    FiniteElementSpace reference_fes;
    DenseTensor local_prolongation;

    real_t h;

    // Builds the local prolongation matrices obtained by refining a single 
    // element.
    void CreateLocalProlongation(int dim);

    // Builds the prolongation matrix from the coarse space to the fine space.
    SparseMatrix CreateProlongation(const VoxelGraph &coarse_graph, const VoxelGraph &fine_graph);
};

// Builds the mesh which contains all neighbor information for a single element.
// In a 2D mesh, this is the 3x3 square. For a 3D mesh, this is the 3x3x3 cube.
Mesh CreateReferenceMesh(int dim);