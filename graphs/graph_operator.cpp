#include "graph_operator.hpp"
#include "mfem.hpp"
#include <iostream>

using namespace mfem;

GraphOperator::GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_,
                             real_t h_)
    : h(h_), reference_fes(reference_fes_), graph(graph_)
{
    // Build A_ref
    BuildARef();

    // Create map from local dof number on central reference element to
    // the local dof number on neigboring reference element
    CreateReferenceDofMap();

    // Create dof_groups and map from broken dofs to true dofs
    CreateDofGroups();

    // Build A matrix
    BuildA(dof_groups);
}

GraphOperator GraphOperator::Coarsen()
{
    Graph coarse_graph = graph.CoarsenGraph();
    return GraphOperator(reference_fes, coarse_graph, h * 2.0);
}

SparseMatrix GraphOperator::CreateProlongation(DenseTensor &local_prolongation,
                                               Graph &fine_graph,
                                               const std::vector<int> &fine_broken_to_true_dof,
                                               const std::vector<std::set<int>> &fine_dof_groups,
                                               const Array<int> &graph_labeling)
{
    static constexpr int map_from_lex[4] = {0, 1, 3, 2};

    const int dofs_per_elem = A_ref.Height();
    const int num_coarse_dofs = dof_groups.size();
    const int num_fine_dofs = fine_dof_groups.size();

    SparseMatrix P(num_fine_dofs, num_coarse_dofs);

    for (int fine_e = 0; fine_e < fine_graph.Size(); ++fine_e)
    {
        const int coarse_e = graph_labeling[fine_e];
        if (coarse_e < 0) { continue; }

        Coord coarse_coord = graph.GetElementCoord(coarse_e);
        Coord fine_coord = fine_graph.GetElementCoord(fine_e);

        // TODO: Generalize to n-dimensions. Currently assumes 2D.
        const int dx = fine_coord[0] - 2 * coarse_coord[0];
        const int dy = fine_coord[1] - 2 * coarse_coord[1];
        if (dx < 0 || dx > 1 || dy < 0 || dy > 1) { MFEM_ABORT(""); }

        const int local_prolongation_index = map_from_lex[dx + 2 * dy];

        Array<int> rows(dofs_per_elem);
        Array<int> cols(dofs_per_elem);
        for (int rf = 0; rf < dofs_per_elem; ++rf)
        {
            const int fine_tdof = fine_broken_to_true_dof[fine_e * dofs_per_elem + rf];
            rows[rf] = fine_tdof;
        }
        for (int rc = 0; rc < dofs_per_elem; ++rc)
        {
            const int coarse_tdof = broken_to_true_dof[coarse_e * dofs_per_elem + rc];
            cols[rc] = coarse_tdof;
        }
        P.SetSubMatrix(rows, cols, local_prolongation(local_prolongation_index));
    }

    // TODO: Handle boundary dofs in P

    P.Finalize();

    return P;
}

void GraphOperator::BuildARef()
{
    DiffusionIntegrator integ1;
    MassIntegrator integ2;

    IsoparametricTransformation T;
    T.SetIdentityTransformation(reference_fes.GetMesh()->GetTypicalElementGeometry());

    auto &PM = T.GetPointMat();
    PM *= h;

    integ1.AssembleElementMatrix(*reference_fes.GetTypicalFE(), T, A_ref);

    DenseMatrix A_ref_2;
    integ2.AssembleElementMatrix(*reference_fes.GetTypicalFE(), T, A_ref_2);
    A_ref += A_ref_2;
}

void GraphOperator::CreateReferenceDofMap()
{
    const int dofs_per_elem = A_ref.Height();
    const int d = reference_fes.GetMesh()->Dimension();

    local_to_neighbor_dof_map.SetSize(dofs_per_elem, pow(3,d));
    local_to_neighbor_dof_map = -1;

    Array<int> ref_dofs;
    int central_element = 0.5 * (pow(3,d)-1);
    reference_fes.GetElementDofs(central_element, ref_dofs);

    // For each local dof of reference element, find the neighboring elements
    // who share the dof, and list the corresponding local dof number on
    // neighbor.
    for (int i = 0; i < dofs_per_elem; ++i)
    {
        int dof = ref_dofs[i];
        for (int e = 0; e < pow(3,d); ++e)
        {
            if (e == central_element) { continue; }
            Array<int> neighbor_dofs;
            reference_fes.GetElementDofs(e, neighbor_dofs);
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
    const int nbdofs = dofs_per_elem * ne;
    Array<int> dof_labeling;

    std::vector<int> element_component_index(nbdofs);
    std::vector<std::set<int>> connected_components(nbdofs);

    // Initialize E_i = { i } and component index
    for (int i = 0; i < nbdofs; ++i)
    {
        connected_components[i].emplace(i);
        element_component_index[i] = i;
    }

    // Method to merge connected component j into connected component i
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
            Array<int> connected_elements = graph.GetConnectedNodes(e);
            for (int neighbor : connected_elements)
            {
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

    for (int i = 0; i < nbdofs; ++i)
    {
        if (!connected_components[i].empty())
        {
            dof_groups.push_back(connected_components[i]);
        }
    }

    // Build map from broken dof number to true dof number
    broken_to_true_dof.assign(nbdofs, -1);
    const int ntdofs = dof_groups.size();
    for (int tdof = 0; tdof < ntdofs; ++tdof)
    {
        const std::set<int> &dof_group = dof_groups[tdof];
        for (int broken_dof : dof_group)
        {
            broken_to_true_dof[broken_dof] = tdof;
        }
    }
}

int GraphOperator::GetNeighborDof(int e, int e_dof, int neighbor)
{
    Coord e_coord = graph.GetElementCoord(e);
    Coord neighbor_coord = graph.GetElementCoord(neighbor);

    int d = reference_fes.GetMesh()->Dimension();
    int neighbor_ref_element = 0;
    for (int i = 0; i < d; ++i)
    {
        neighbor_ref_element += (neighbor_coord[i] - e_coord[i] + 1) * pow(3, i);
    }
    return local_to_neighbor_dof_map(e_dof, neighbor_ref_element);
}

void GraphOperator::BuildA(std::vector<std::set<int>> dof_groups)
{
    // This method feels very unoptimized. The matrices are very sparse and
    // have a lot of structure, which could be explioted. However, this should
    // work, and was the most direct way I could see using MFEM's methods.

    // Build A using matrix. Lambda(i,j) = 1 if broken dof i is identified with
    // true dof j, and 0 otherwise.

    const int dofs_per_elem = A_ref.Height();
    const int ne = graph.Size();
    const int nbdofs = dofs_per_elem * ne;
    const int ntdofs = dof_groups.size();

    SparseMatrix Lambda(nbdofs, ntdofs);
    for (int bdof = 0; bdof < nbdofs; ++bdof)
    {
        int tdof = broken_to_true_dof[bdof];
        Lambda.Set(bdof, tdof, 1.0);
    }

    // Build hat{A} matrix.
    SparseMatrix A_hat(nbdofs, nbdofs);

    for (int e = 0; e < ne; ++e)
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

    {
        std::ofstream f("A_hat.txt");
        A_hat.PrintMatlab(f);
    }

    // Compute A = Lambda^T * hat{A} * Lambda
    A = *RAP(Lambda, A_hat, Lambda);
    A.Finalize();

    // TODO: Eliminate boundary dofs from coarse operator
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
