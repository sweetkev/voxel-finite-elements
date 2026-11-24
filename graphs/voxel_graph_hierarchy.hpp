# pragma once

#include "mfem.hpp"
#include "voxel_graph.hpp"
#include "voxel_graph_operator.hpp"
#include "cell_ownership.hpp"

using namespace mfem;
using namespace std;

class VoxelGraphHierarchy
{
public:
    // The finest element matrices are built from a's integrators,
    // and each coarser level operator from the element matrices of the next
    // finer level. ownership gives the rank that coarsens each coarse grid
    // cell.
    VoxelGraphHierarchy(unique_ptr<VoxelGraph> fine_graph, Array<int> &fine_ess_dofs, int nlevels,
                        const CellOwnership &ownership,
                        FiniteElementSpace &reference_fes, real_t h,
                        ParBilinearForm &a, HypreParMatrix &fine_A);

    vector<unique_ptr<VoxelGraphOperator>> &GetGraphOperators() { return graph_operators; }

    // Parallel operator on the given level (0 is coarsest). On the finest
    // level this is fine_A.
    HypreParMatrix &GetParOperator(int level)
    {
        if (level == (int)graph_operators.size() - 1) { return fine_A; }
        return graph_operators[level]->GetParMatrix();
    }

    // Parallel prolongation from level to level + 1 (true dofs to true dofs),
    // with the essential true dofs of both levels eliminated.
    HypreParMatrix &GetParProlongation(int level)
    {
        return *prolongations[level];
    }

    const VoxelGraph &GetGraph(int level) const { return *graphs[level]; }

    // Essential true dofs (rank-local true dof indices) on the given level.
    const Array<int> &GetEssentialTrueDofs(int level) const
    {
        return ess_true_dofs[level];
    }

private:
    vector<unique_ptr<VoxelGraph>> graphs;
    // Essential true dofs (rank-local true dof indices) on each level
    vector<Array<int>> ess_true_dofs;
    // Prolongations from level i to level i + 1 (true dof numbering)
    vector<unique_ptr<HypreParMatrix>> prolongations;
    vector<unique_ptr<VoxelGraphOperator>> graph_operators;
    // Finest level operator, assembled by MFEM
    HypreParMatrix &fine_A;

    // Owner rank of each coarse grid cell on every coarse level
    const CellOwnership &ownership;

    FiniteElementSpace &reference_fes;
    DenseTensor local_prolongation;

    real_t h;

    // Builds the local prolongation matrices obtained by refining a single
    // element.
    void CreateLocalProlongation(int dim);

    // Builds the prolongation from the coarse true dofs to the fine true dofs
    // as P = D^{-1} * Lambda_f^T * P-hat * Lambda_c, where Lambda_c and
    // Lambda_f are the graphs' broken-to-true-dof matrices, P-hat is the
    // broken prolongation (one local prolongation block per fine element), and
    // D counts the broken copies of each fine true dof. Sets
    // coarse_ess_true_dofs from fine_ess_true_dofs; the essential rows and
    // columns are not eliminated from P. Collective.
    unique_ptr<HypreParMatrix> CreateProlongation(
        const VoxelGraph &coarse_graph, const VoxelGraph &fine_graph,
        Array<int> &coarse_ess_true_dofs, const Array<int> &fine_ess_true_dofs);
};

// Builds the mesh which contains all neighbor information for a single element.
// In a 2D mesh, this is the 3x3 square. For a 3D mesh, this is the 3x3x3 cube.
Mesh CreateReferenceMesh(int dim);
