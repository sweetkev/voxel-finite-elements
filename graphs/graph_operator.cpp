#include "graph_operator.hpp"
#include "mfem.hpp"

using namespace mfem;

GraphOperator::GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_)
    : reference_fes(reference_fes_), graph(graph_)
{
    // Create map from local dof # to neighboring cells with dof
    CreateDofToElementMap();

    // Create global dof labeling
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
        for (int i = 0; i < dofs_per_elem; ++i)
        {
            Array<int> neighbors_with_dof = NeighborsWithDof(e,i);
        }
    }
}

void GraphOperator::CreateDofToElementMap()
{
    Array<int> ref_dofs;
    reference_fes.GetElementDofs(4, ref_dofs);
    const int dofs_per_elem = ref_dofs.Size();

    Table dof_to_element_map;

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
                    dof_to_element_map.Push(i, e);
                    break;
                }
            }
        }
    }

    dof_to_reference_element_map = dof_to_element_map;
}

Array<int> GraphOperator::NeighborsWithDof(int e, int dof)
{
    Array<int> neighbors_with_dof;
}

void CreateReferenceMesh(Mesh &reference_mesh, int dim) 
{
    /** Creates 3x3 or 3x3x3 reference mesh */
    //TODO: implement 3-dimensions
    //MFEM_ASSERT(dim == 2 || dim == 3, "Dimension must be 2 or 3");
        
    MFEM_ASSERT(dim == 2, "dimension must be 2");

    Mesh mesh;
    /** 3x3 reference mesh with elements and vertices numbered as follows:
     * elements:    vertices:
     *              12 13 14 15
     *  6 7 8       8  9  10 11
     *  3 4 5       4  5  6  7
     *  0 1 2       0  1  2  3
     */
    if(dim == 2)
    {
        for(int i = 0; i < 4; i++) 
        {
            for(int j = 0; j < 4; j++) 
            {
                mesh.AddVertex(i,j);
            }
        }

        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++) 
            {
                mesh.AddQuad(i+j*4, i+1+j*4, i+1+(j+1)*4, i+(j+1)*4);
            }
        }
    }

    mesh.Finalize();
    
    Swap(mesh, reference_mesh);
}