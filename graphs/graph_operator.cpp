#include "graph_operator.hpp"
#include "mfem.hpp"

using namespace mfem;

GraphOperator::GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_)
    : reference_fes(reference_fes_), graph(graph_)
{
    // Create global dof labeling
    Array<int> ref_dofs;
    reference_fes.GetElementDofs(0, ref_dofs);
    int dofs_per_elem = ref_dofs.Size();
    int ne = graph.Size();

    Array<int> dof_labeling;
    dof_labeling.SetSize(ne * dofs_per_elem, -1);
    int dof_idx = 0;
    for(int elem = 0; elem < ne; ++elem)
    {
        for(int i = 0; i < dofs_per_elem; ++i)
        {
            if(dof_labeling[elem * dofs_per_elem + i] == -1)
            {
                dof_labeling[elem * dofs_per_elem + i] = dof_idx;
                LabelSharedDofs(elem, i, dof_idx, dof_labeling, dofs_per_elem);
                dof_idx++;    
            }
        }
    }
    int ndofs = dof_idx;

}

void GraphOperator::LabelSharedDofs(int elem, int dof, int dof_idx,
                         Array<int> &dof_labeling, int dofs_per_elem)
{ 
    Array<int> elem_dofs;
    if(reference_fes.GetMesh()->Dimension() == 2)
    {
        // elements in reference mesh that contain the given dof
        Array<int> reference_elements;
        reference_elements.SetSize(9, -1);
        // reference_elements_dofs[i] contains the local dof number corresponding to
        // the dof with respect to reference_elements[i]
        Array<int> reference_elements_dofs;
        reference_elements_dofs.SetSize(9, -1);

        reference_fes.GetElementDofs(4, elem_dofs);
        int global_dof = elem_dofs[dof];
        for(int i = 0; i < 9; ++i)
        {
            if(i != 4)
            {
                reference_fes.GetElementDofs(i, elem_dofs);
                for(int j = 0; j < elem_dofs.Size(); ++j)
                {
                    if(elem_dofs[j] == global_dof)
                    {
                        reference_elements.Append(i);
                        reference_elements_dofs.Append(j);
                        break;
                    }
                }
            }
        }

        // If the dof is an interior dof, then it is not shared and we can 
        // return early
        bool interior_dof = true;
        for(int i = 0; i < reference_elements.Size(); ++i)
        {
            if(reference_elements[i] != -1)
            {
                interior_dof = false;
                break;
            }
        }
        if(interior_dof) { return; }
    
        // List of elements that share given dof
        Array<int> elements_with_dof;
        // CSR-style pointer array to store the element type of the elements 
        // that share the given dof
        Array<int> element_pointers;
        element_pointers.SetSize(10,0);

        GetElementsWithDof(elem, dof, reference_elements, 
            reference_elements_dofs, elements_with_dof, element_pointers);

        // Label all local dofs that are identified with the current dof
        for(int i = 0; i < element_pointers.Size(); ++i)
        {
            for(int j = element_pointers[i]; j < element_pointers[i+1]; ++j)
            {
                int k = elements_with_dof[j] * dofs_per_elem 
                        + reference_elements_dofs[i];
                dof_labeling[k] = dof_idx;
            }   
        }
    }
}

void GraphOperator::GetElementsWithDof(int elem, int dof, 
    Array<int> &reference_elements, Array<int> &reference_elements_dofs, 
    Array<int> &elements_with_dof, Array<int> &element_pointers)
{
    if(reference_fes.GetMesh()->Dimension() == 2)
    {
        // element_groups[i] contains the list of elements that that share the
        // dof with current element, and correspond to element i in the
        // reference mesh
        std::vector<std::vector<int>> element_groups(9);

        Array<int> connected_nodes = graph.GetConnectedNodes(elem);
        for(int i = 0; i < connected_nodes.Size(); ++i)
        {
            int connected_elem = connected_nodes[i];
            int connected_cell = graph.GetNodeCell(connected_elem);
            Coord connected_cell_coord = graph.GetCellCoord(connected_cell);
            int elem_cell = graph.GetNodeCell(elem);
            Coord elem_cell_coord = graph.GetCellCoord(elem_cell);

            // Find which element in the reference mesh corresponds to the
            // connected element, and if it shares the dof with current element
            switch (connected_cell_coord[0] - elem_cell_coord[0])
            {
                // Element to left
                case -1:
                    switch (connected_cell_coord[1] - elem_cell_coord[1])
                    {
                        // Bottom-left element
                        case -1:
                            if(reference_elements[0] >= 0)
                            {
                                element_groups[0].push_back(connected_elem);
                            }
                            break;
                        // Middle-left element
                        case 0:
                            if(reference_elements[1] >= 0)
                            {
                                element_groups[1].push_back(connected_elem);
                            }
                            break;
                        // Top-left element
                        case 1:
                            if(reference_elements[2] >= 0)
                            {
                                element_groups[2].push_back(connected_elem);
                            }
                            break;
                        default:
                            break;
                    }
                    break;
                // Element in same column
                case 0:
                    switch (connected_cell_coord[1] - elem_cell_coord[1])
                    {
                        // Bottom-middle element
                        case -1:
                            if(reference_elements[3] >= 0)
                            {
                                element_groups[3].push_back(connected_elem);
                            }
                            break;
                        // Middle element
                        case 0:
                            if(reference_elements[4] >= 0)
                            {
                                element_groups[4].push_back(connected_elem);
                            }
                            break;
                        // Top-middle element
                        case 1:
                            if(reference_elements[5] >= 0)
                            {
                                element_groups[5].push_back(connected_elem);
                            }
                            break;
                        default:
                            break;
                    }
                    break;
                // Element to right
                case 1:
                    switch (connected_cell_coord[1] - elem_cell_coord[1])
                    {
                        // Bottom-right element
                        case -1:
                            if(reference_elements[6] >= 0)
                            {
                                element_groups[6].push_back(connected_elem);
                            }
                            break;
                        // Middle-right element
                        case 0:
                            if(reference_elements[7] >= 0)
                            {
                                element_groups[7].push_back(connected_elem);
                            }
                            break;
                        // Top-right element
                        case 1:
                            if(reference_elements[8] >= 0)
                            {
                                element_groups[8].push_back(connected_elem);
                            }
                            break;
                        default:
                            break;
                    }
                    break;
                default:
                    break;
            }

            // Place all elements that share the dof in a single list, and
            // add pointers to the different element groups
            for(int j = 0; j < reference_elements.Size(); ++j)
            {
                for(int k = 0; k < element_groups[j].size(); ++k)
                {
                    elements_with_dof.Append(element_groups[j][k]);
                }
                element_pointers[j+1] = element_pointers[j] + 
                                        element_groups[j].size();
            }
        }
    }    
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