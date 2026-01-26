#include "graph_operator.hpp"
#include "mfem.hpp"

using namespace mfem;

void GraphOperator::CreateReferenceMesh(int dim) 
{
    /** Creates 3x3 or 3x3x3 reference mesh */
    //TODO: implement 3-dimensions
    //MFEM_ASSERT(dim == 2 || dim == 3, "Dimension must be 2 or 3");
        
    MFEM_ASSERT(dim == 2, "dimension must be 2");

    Mesh mesh;

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
            for(int j = 0; j < 3; i++) 
            {
                mesh.AddQuad(i+j*4, i+1+j*4, i+1+(j+1)*4, i+(j+1)*4);
            }
        }
    }

    mesh.Finalize();
    reference_mesh.Swap(mesh,true); 
}