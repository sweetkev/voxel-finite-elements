#include "graph_operator.hpp"
#include "mfem.hpp"

using namespace mfem;

GraphOperator::GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_)
    : reference_fes(reference_fes_), graph(graph_)
{
    // Create map from local dof # to neighboring cells with dof
    CreateDofMap();

    // Create dof_groups
    std::vector<std::set<int>> dof_groups = CreateDofGroups();

    // Build Q matrix
    BuildQ(dof_groups);
}

void GraphOperator::CreateDofMap()
{
    Array<int> ref_dofs;
    reference_fes.GetElementDofs(4, ref_dofs);
    const int dofs_per_elem = ref_dofs.Size();

    local_to_neighbor_dof_map.SetSize(dofs_per_elem, 9);

    // Initialize local dof - to - neighbor dof matrix with -1
    for (int i = 0; i < local_to_neighbor_dof_map.Width(); ++i)
    {
        for (int j = 0; j < local_to_neighbor_dof_map.Height(); ++j)
        {
            local_to_neighbor_dof_map(i,j) = -1;
        }
    }

    // For each DoF of reference element, find the neighboring elements who
    // share the Dof
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
    Array<int> ref_dofs;
    reference_fes.GetElementDofs(0, ref_dofs);
    const int dofs_per_elem = ref_dofs.Size();
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
            switch (neighbor_coord[0] - neighbor_coord[0])
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
                    break;
            }
        // Middle elements
        case 0:
            switch (neighbor_coord[1] - neighbor_coord[1])
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
                    break;
            }
        // Top elements
        case 1:
            switch (neighbor_coord[1] - neighbor_coord[1])
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
                    break;
            }

        default:
            return -1;
    }
}

void GraphOperator::BuildQ(std::vector<std::set<int>> dof_groups)
{

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
