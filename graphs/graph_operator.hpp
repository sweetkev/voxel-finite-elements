#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include "graph.hpp"
#include <unordered_map>

using namespace mfem;

class GraphOperator
{
public:
    GraphOperator(FiniteElementSpace &reference_fes_, Graph &graph_);

private:
    Array<int> A_ref;
    FiniteElementSpace reference_fes;
    SparseMatrix Q;
    Graph graph;

    /** Assigns the label 'dof_idx' to all local DoFs that are identified with
     * the current local DoF of element 'elem'.
     */
    void LabelSharedDofs(int elem, int local_dof, int dof_idx,
                         Array<int> &dof_labeling, int dofs_per_elem);

    /** Returns the list of nodes that share the given local Dof of element
     * 'elem'.
     */
    void GraphOperator::GetElementsWithDof(int elem, int dof, 
        Array<int> &reference_elements, Array<int> &reference_elements_dofs, 
        Array<int> &elements_with_dof, Array<int> &element_pointers);
};

void CreateReferenceMesh(Mesh &reference_mesh,int dim);