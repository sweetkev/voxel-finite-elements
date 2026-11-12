# pragma once

#include "mfem.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"

using namespace mfem;
using namespace std;

class VoxelGraphHierarchy
{
public:
    // Creates finest level operator from a.
    // TODO: Implement coarser levels.
    VoxelGraphHierarchy(unique_ptr<VoxelGraph> fine_graph, Array<int> &fine_ess_dofs, int nlevels,
                        FiniteElementSpace &reference_fes, real_t h,
                        ParBilinearForm &a);

    vector<unique_ptr<VoxelGraphOperator>> &GetGraphOperators() { return graph_operators; }
    vector<unique_ptr<SparseMatrix>> &GetProlongations() { return prolongations; }

    // Parallel operator on the given level (0 is coarsest).
    HypreParMatrix &GetParOperator(int level)
    {
        return graph_operators[level]->GetParMatrix();
    }

    // Parallel prolongation from level to level + 1 (true dofs to true dofs).
    HypreParMatrix &GetParProlongation(int level)
    {
        return *par_prolongations[level];
    }

private:
    vector<unique_ptr<VoxelGraph>> graphs;
    vector<Array<int>> ess_dofs;
    // Rank-local prolongations (local dof numbering)
    vector<unique_ptr<SparseMatrix>> prolongations;
    // Parallel prolongations (true dof numbering)
    vector<unique_ptr<HypreParMatrix>> par_prolongations;
    vector<unique_ptr<VoxelGraphOperator>> graph_operators;
    

    FiniteElementSpace &reference_fes;
    DenseTensor local_prolongation;

    real_t h;

    // Builds the local prolongation matrices obtained by refining a single
    // element.
    void CreateLocalProlongation(int dim);

    // Builds the prolongation matrix from the coarse space to the fine space.
    SparseMatrix CreateProlongation(const VoxelGraph &coarse_graph, 
                                    const VoxelGraph &fine_graph, 
                                    Array<int> &coarse_ess_dofs,
                                    const Array<int> &fine_ess_dofs);
};

// Builds the mesh which contains all neighbor information for a single element.
// In a 2D mesh, this is the 3x3 square. For a 3D mesh, this is the 3x3x3 cube.
Mesh CreateReferenceMesh(int dim);
