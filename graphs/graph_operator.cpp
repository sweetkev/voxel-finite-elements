#include "graph_operator.hpp"
#include "mfem.hpp"

using namespace mfem;

GraphOperator::GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_)
    : reference_fes(reference_fes_), graph(graph_)
{
    // Build A_ref
    BuildARef();

    // Create map from local dof # to neighboring cells with dof
    CreateDofMap();

    // Create dof_groups
    std::vector<std::set<int>> dof_groups = CreateDofGroups();

    // Build A matrix
    BuildA(dof_groups);
}

GraphOperator GraphOperator::Coarsen()
{
    Graph coarse_graph = graph.CoarsenGraph();
    return GraphOperator(A_ref, reference_fes, coarse_graph);
}

GraphOperator::GraphOperator(DenseMatrix &A_ref_, FiniteElementSpace &reference_fes_, 
                  Graph &graph_)
    : A_ref(A_ref_), reference_fes(reference_fes_), graph(graph_)
{
    // Create map from local dof # to neighboring cells with dof
    CreateDofMap();

    // Create dof_groups
    std::vector<std::set<int>> dof_groups = CreateDofGroups();

    // Build A matrix
    BuildA(dof_groups);
}

SparseMatrix GraphOperator::CreateProlongation(DenseTensor &local_prolongation,
                                               Graph &fine_graph)
{
    int dofs_per_elem = A_ref.Height();
    SparseMatrix P(fine_graph.Size()*dofs_per_elem,
                   graph.Size()*dofs_per_elem);

    for (int i = 0; i < graph.Size(); ++i)
    {
        // TODO: Generalize to n-dimensions
        // Currently assumes 2D
        Coord e_coord = graph.GetElementCoord(i);
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                Coord fine_coord(e_coord[0] + j, e_coord[1] + k);
                const std::vector<Coord> &fine_grid_cells = fine_graph.GetGridCells();
                auto it = std::find(fine_grid_cells.begin(), fine_grid_cells.end(), 
                                    fine_coord);
            
                if (it != fine_grid_cells.end())
                {
                    int neighbor_index = std::distance(fine_grid_cells.begin(), it);

                    // prepare row/column index arrays for the submatrix.
                    // rows correspond to fine-graph dofs, cols to coarse-graph dofs.
                    const int nd = local_prolongation.SizeI();
                    const int nc = local_prolongation.SizeJ();
                    Array<int> rows(nd);
                    Array<int> cols(nc);

                    for (int r = 0; r < nd; ++r)
                    {
                        rows[r] = nd * neighbor_index + r;
                    }

                    for (int c = 0; c < nc; ++c)
                    {
                        cols[c] = nc * i + c;
                    }

                    P.AddSubMatrix(rows, cols, local_prolongation(j + 2*k));
                }
            }
        }    
    }

    RemoveBoundaryDofs(P);

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

std::vector<std::set<int>> GraphOperator::CreateDofGroups()
{
    const int dofs_per_elem = A_ref.Height();
    const int ne = graph.Size();
    const int nd = dofs_per_elem * ne;
    Array<int> dof_labeling;

    std::vector<int> element_component_index(nd);
    std::vector<std::set<int>> connected_components(nd);

    // Initialize E_i = { i }
    for (int i = 0; i < nd; ++i)
    {
        connected_components[i].emplace(i);
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

    std::vector<std::set<int>> dof_groups;
    for (int i = 0; i < nd; ++i)
    {
        if (!connected_components[i].empty())
        {
            dof_groups.push_back(connected_components[i]);
        }
    }

    return dof_groups;
}

int GraphOperator::GetNeighborDof(int e, int e_dof, int neighbor)
{
    Coord e_coord = graph.GetElementCoord(e);
    Coord neighbor_coord = graph.GetElementCoord(neighbor);

    switch (neighbor_coord[1] - e_coord[1])
    {
        // Bottom elements
        case -1:
            switch (neighbor_coord[0] - e_coord[0])
            {
                // Bottom-left element
                case -1:
                    return local_to_neighbor_dof_map(e_dof, 0);
                // Bottom-middle element
                case 0:
                    return local_to_neighbor_dof_map(e_dof, 1);
                // Bottom-right element
                case 1:
                    return local_to_neighbor_dof_map(e_dof, 2);
                default:
                    return -1;
            }
        // Middle elements
        case 0:
            switch (neighbor_coord[0] - e_coord[0])
            {
                // Middle-left element
                case -1:
                    return local_to_neighbor_dof_map(e_dof, 3);
                // Middle-middle element
                case 0:
                    return local_to_neighbor_dof_map(e_dof, 4);
                // Middle-right element
                case 1:
                    return local_to_neighbor_dof_map(e_dof, 5);
                default:
                    return -1;
            }
        // Top elements
        case 1:
            switch (neighbor_coord[0] - e_coord[0])
            {
                // Top-left element
                case -1:
                    return local_to_neighbor_dof_map(e_dof, 6);
                // Top-middle element
                case 0:
                    return local_to_neighbor_dof_map(e_dof, 7);
                // Top-right element
                case 1:
                    return local_to_neighbor_dof_map(e_dof, 8);
                default:
                    return -1;
            }

        default:
            return -1;
    }
}

void GraphOperator::BuildA(std::vector<std::set<int>> dof_groups)
{
    // This method feels very unoptimized. The matrices are very sparse and
    // have a lot of structure, which could be explioted. However, this should
    // work, and was the most direct way I could see using MFEM's existing
    // methods.

    // Build A using matrix. Lambda(i,j) = 1 if broken dof i is identified with
    // global dof j, and 0 otherwise.
    
    const int dofs_per_elem = A_ref.Height();
    const int ne = graph.Size();
    const int total_dofs = dofs_per_elem * ne;
    
    SparseMatrix Lambda(total_dofs, dof_groups.size());
    for (int global_dof = 0; global_dof < dof_groups.size(); ++global_dof)
    {
        std::set<int> group = dof_groups[global_dof];
        for(int broken_dof : group)
        {
            Lambda.Add(broken_dof, global_dof, 1.0);
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
    A = *RAP(A_hat, Lambda);
}

void GraphOperator::RemoveBoundaryDofs(SparseMatrix &P)
{ 
    const int dofs_per_elem = A_ref.Height();
    for (int e = 0; e < graph.Size(); ++e)
    {
        Array<int> connected_elements = graph.GetConnectedNodes(e);
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

            // Check if dof is on boundary
            Coord e_coord = graph.GetElementCoord(e);
            Array<int> neighbors_to_check(9);
            neighbors_to_check = -1;
            for (int i = 0; i < 9; ++i)
            {
                if (local_to_neighbor_dof_map(e_dof,i) != -1)
                {
                    neighbors_to_check[i] = 1;
                }
            }

            for (int element : connected_elements)
            {
                Coord neighbor_coord = graph.GetElementCoord(element);
                int neighbor_index = (neighbor_coord[0] - e_coord[0] + 1)
                                     + 3*(neighbor_coord[1] - e_coord[1] + 1);
                neighbors_to_check[neighbor_index] = -1;
            }

            for (int neighbor : neighbors_to_check)
            {
                if (neighbor == 1)
                {
                    P.EliminateRow(e * dofs_per_elem + e_dof);

                    // Also remove boundary condition from A. This should
                    // probably be implemented elsewhere, but this is most
                    // efficient for now.
                    A.EliminateRow(e * dofs_per_elem + e_dof, Operator::DIAG_ONE);
                }
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
