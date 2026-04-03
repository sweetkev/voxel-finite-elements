#include "graph_operator.hpp"
#include "mfem.hpp"
#include <iostream>

using namespace mfem;

GraphOperator::GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_)
    : reference_fes(reference_fes_), graph(graph_)
{
    // Build A_ref
    BuildARef();

    // Create map from local dof # to neighboring cells with dof
    CreateDofMap();

    // Create dof_groups
    CreateDofGroups();

    // Build A matrix
    BuildA(dof_groups);
}

GraphOperator GraphOperator::Coarsen()
{
    Graph coarse_graph = graph.CoarsenGraph();
    return GraphOperator(A_ref, reference_fes, coarse_graph, level + 1);
}

GraphOperator::GraphOperator(DenseMatrix A_ref_, FiniteElementSpace &reference_fes_, 
                  Graph &graph_, int level_)
    : A_ref(A_ref_), reference_fes(reference_fes_), graph(graph_), level(level_)
{
    A_ref.Set(1 / pow(pow(2, level), reference_fes.GetFE(0)->GetDim()), A_ref);
    // Create map from local dof # to neighboring cells with dof
    CreateDofMap();

    // Create dof_groups
    CreateDofGroups();

    // Build A matrix
    BuildA(dof_groups);
}

SparseMatrix GraphOperator::CreateProlongation(DenseTensor &local_prolongation,
                                               Graph &fine_graph,
                                               const std::vector<int> &fine_broken_to_true_dof,
                                               const std::vector<std::set<int>> fine_dof_groups)
{
    const int dofs_per_elem = A_ref.Height();
    const int num_coarse_dofs = dof_groups.size();
    const int num_fine_dofs = fine_dof_groups.size();

    SparseMatrix P(num_fine_dofs, num_coarse_dofs);

    // Build prolongation in true DOF space using the broken->true mapping.
    for (int i = 0; i < graph.Size(); ++i)
    {
        // TODO: Generalize to n-dimensions
        // Currently assumes 2D
        Coord e_coord = graph.GetElementCoord(i);
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                Coord fine_coord(2 * e_coord[0] + j, 2 * e_coord[1] + k);
                const std::vector<Coord> &fine_grid_cells = fine_graph.GetGridCells();
                auto it = std::find(fine_grid_cells.begin(), fine_grid_cells.end(), 
                                    fine_coord);
            
                if (it != fine_grid_cells.end())
                {
                    int neighbor_index = std::distance(fine_grid_cells.begin(), it);

                    int coarse_base = i * dofs_per_elem;
                    int fine_base = neighbor_index * dofs_per_elem;

                    // Insert into P using true global dof indices.
                    Array<int> rows(dofs_per_elem);
                    Array<int> cols(dofs_per_elem);
                    for (int rf = 0; rf < dofs_per_elem; ++rf)
                    {
                        const int fine_tdof = fine_broken_to_true_dof[fine_base + rf];
                        rows[rf] = fine_tdof;
                    }
                    for (int rc = 0; rc < dofs_per_elem; ++rc)
                    {
                        const int coarse_tdof = broken_dof_to_true_dof[coarse_base + rc];
                        cols[rc] = coarse_tdof;
                    }
                    P.AddSubMatrix(rows, cols, local_prolongation(j + 2*k));
                }
            }
        }
    }

    // Remove boundary dofs (now correctly implemented in true DOF space)
    RemoveBoundaryDofs(P, Operator::DIAG_ZERO);

    P.Finalize();

    return P;
}

void GraphOperator::BuildARef()
{
    Array<int> *ess_dofs = new Array<int>;
    reference_fes.GetBoundaryTrueDofs(*ess_dofs);
    BilinearForm a(&reference_fes);
    a.AddDomainIntegrator(new DiffusionIntegrator);
    a.Assemble();
    SparseMatrix large_A = a.SpMat();

    Array<int> ref_dofs;
    reference_fes.GetElementDofs(4, ref_dofs);
    
    const int dofs_per_elem = ref_dofs.Size();
    A_ref.SetSize(dofs_per_elem, dofs_per_elem);
    
    large_A.GetSubMatrix(ref_dofs, ref_dofs, A_ref);
}

void GraphOperator::CreateDofMap()
{
    const int dofs_per_elem = A_ref.Height();

    local_to_neighbor_dof_map.SetSize(dofs_per_elem, 9);
    local_to_neighbor_dof_map = -1;

    // For each DoF of reference element, find the neighboring elements who
    // share the Dof
    Array<int> ref_dofs;
    reference_fes.GetElementDofs(4, ref_dofs);
    for (int i = 0; i < dofs_per_elem; ++i)
    {
        int dof = ref_dofs[i];
        for (int e = 0; e < 9; ++e)
        {
            if (e == 4) { continue; }
            Array<int> neighbor_dofs;
            reference_fes.GetElementDofs(e,neighbor_dofs);
            for (int j = 0; j < neighbor_dofs.Size(); ++j)
            {
                if (neighbor_dofs[j] == dof)
                {
                    local_to_neighbor_dof_map(i,e) = j;
                    break;
                }
            }
        }
    }
}

void GraphOperator::CreateDofGroups()
{
    const int dofs_per_elem = A_ref.Height();
    const int ne = graph.Size();
    const int nd = dofs_per_elem * ne;
    Array<int> dof_labeling;

    std::vector<int> element_component_index(nd);
    std::vector<std::set<int>> connected_components(nd);

    // Initialize E_i = { i } and component index
    for (int i = 0; i < nd; ++i)
    {
        connected_components[i].emplace(i);
        element_component_index[i] = i;
    }

    // merges connected component j into connected component i
    auto merge = [&](int i, int j)
    {
        const int ci = element_component_index[i];
        const int cj = element_component_index[j];
        element_component_index[j] = ci;
        connected_components[ci].merge(connected_components[cj]);
    };

    for (int e = 0; e < ne; ++e)
    {
        for (int e_dof = 0; e_dof < dofs_per_elem; ++e_dof)
        {
            // Skip if interior dof
            bool interior_dof = true;
            for (int j = 0; j < 9; ++j)
            {
                if (local_to_neighbor_dof_map(e_dof,j) != -1)
                {
                    interior_dof = false;
                }
            }
            if (interior_dof) { continue; }

            // For each connected element, check if the dof is shared between
            // them.
            Array<int> connected_elements = graph.GetConnectedNodes(e);
            for (int j = 0; j < connected_elements.Size(); ++j)
            {
                int neighbor = connected_elements[j];
                int neighbor_dof = GetNeighborDof(e, e_dof, neighbor);
                if (neighbor_dof != -1)
                {
                    int e_index = dofs_per_elem * e + e_dof;
                    int neighbor_index = dofs_per_elem * neighbor 
                                         + neighbor_dof;
                    int ci = element_component_index[e_index];
                    int cj = element_component_index[neighbor_index];
                    if (ci != cj)
                    {
                        merge(e_index, neighbor_index);
                    }
                }
            }
        }
    }

    for (int i = 0; i < nd; ++i)
    {
        if (!connected_components[i].empty())
        {
            dof_groups.push_back(connected_components[i]);
        }
    }
}

int GraphOperator::GetNeighborDof(int e, int e_dof, int neighbor)
{
    Coord e_coord = graph.GetElementCoord(e);
    Coord neighbor_coord = graph.GetElementCoord(neighbor);

    int neighbor_ref_element = (neighbor_coord[0] - e_coord[0] + 1) 
                             + 3 * (neighbor_coord[1] - e_coord[1] + 1);
    return local_to_neighbor_dof_map(e_dof, neighbor_ref_element);
}

void GraphOperator::BuildA(std::vector<std::set<int>> dof_groups)
{
    // This method feels very unoptimized. The matrices are very sparse and
    // have a lot of structure, which could be explioted. However, this should
    // work, and was the most direct way I could see using MFEM's existing
    // methods.

    // Build A using matrix. Lambda(i,j) = 1 if broken dof i is identified with
    // true global dof j, and 0 otherwise.
    
    const int dofs_per_elem = A_ref.Height();
    const int ne = graph.Size();
    const int total_dofs = dofs_per_elem * ne;

    // Map from broken dof to global dof (true dof numbering)
    broken_dof_to_true_dof.assign(total_dofs, -1);

    SparseMatrix Lambda(total_dofs, dof_groups.size());
    for (int global_dof = 0; global_dof < dof_groups.size(); ++global_dof)
    {
        const std::set<int> &group = dof_groups[global_dof];
        for (int broken_dof : group)
        {
            Lambda.Add(broken_dof, global_dof, 1.0);
            broken_dof_to_true_dof[broken_dof] = global_dof;
        }
    }

    // Build hat{A} matrix. (This feels like a very silly way of doing this,
    // but it should work)
    SparseMatrix A_hat(total_dofs, total_dofs);
    
    for (int e = 0; e < graph.Size(); ++e)
    {
        Array<int> block_indices(dofs_per_elem);
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            block_indices[i] = e * dofs_per_elem + i;
        }
        A_hat.AddSubMatrix(block_indices, block_indices, A_ref);
    }

    // Finalize matrices before RAP operation
    Lambda.Finalize();
    A_hat.Finalize();

    // Compute A = Lambda^T * hat{A} * Lambda
    A = *RAP(Lambda, A_hat, Lambda);

    // Eliminate boundary dofs from coarse operator
    RemoveBoundaryDofs(A, Operator::DIAG_ONE);
}

void GraphOperator::RemoveBoundaryDofs(SparseMatrix &B, 
                                       Operator::DiagonalPolicy dpolicy)
{ 
    const int dofs_per_elem = A_ref.Height();
    const int ne = graph.Size();

    // Identify and eliminate boundary dofs
    for (int e = 0; e < ne; ++e)
    {
        Array<int> connected_elements = graph.GetConnectedNodes(e);
        
        for (int e_dof = 0; e_dof < dofs_per_elem; ++e_dof)
        {
            // Skip if interior dof
            bool interior_dof = true;
            for (int j = 0; j < 9; ++j)
            {
                if (local_to_neighbor_dof_map(e_dof, j) != -1)
                {
                    interior_dof = false;
                    break;
                }
            }
            if (interior_dof) { continue; }

            // Check if this dof is on the boundary of the domain
            Coord e_coord = graph.GetElementCoord(e);
            bool on_boundary = false;

            // Check all 8 potential neighbor directions (skip center element 4)
            for (int neighbor_dir = 0; neighbor_dir < 9; ++neighbor_dir)
            {
                if (neighbor_dir == 4) continue;
                
                if (local_to_neighbor_dof_map(e_dof, neighbor_dir) != -1)
                {
                    // Convert direction index to coordinate offset
                    int dx = (neighbor_dir % 3) - 1;  // -1, 0, or 1
                    int dy = (neighbor_dir / 3) - 1;  // -1, 0, or 1
                    Coord neighbor_coord(e_coord[0] + dx, e_coord[1] + dy);

                    // Check if this neighbor exists in the graph
                    bool neighbor_exists = false;
                    for (int neighbor : connected_elements)
                    {
                        if (graph.GetElementCoord(neighbor) == neighbor_coord)
                        {
                            neighbor_exists = true;
                            break;
                        }
                    }

                    if (!neighbor_exists)
                    {
                        on_boundary = true;
                        break;
                    }
                }
            }

            if (on_boundary)
            {
                int broken_dof_index = e * dofs_per_elem + e_dof;
                int true_dof_index = broken_dof_to_true_dof[broken_dof_index];
                
                B.EliminateRow(true_dof_index, dpolicy);
            }
        }
    }
}

Mesh CreateReferenceMesh(int dim)
{
    if (dim == 2)
    {
        return Mesh::MakeCartesian2D(3, 3, Element::QUADRILATERAL, false, 
                                     1.0, 1.0, false);
    }
    else if (dim == 3)
    {
        return Mesh::MakeCartesian3D(3, 3, 3, Element::HEXAHEDRON, 
                                     1.0, 1.0, 1.0, false);
    }
    else
    {
        MFEM_ABORT("Unsupported dimension.");
    }
}
