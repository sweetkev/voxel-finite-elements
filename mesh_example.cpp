//ensure path to mfem.hpp is correct
#include "/home/kevin/mfem/mfem-4.8/mfem.hpp"

using namespace mfem;

int main() {
    int dim = 2, nv = 9, ne = 4, nb = 0, sdim = 2, order = 1, attrib = 1;
    Mesh mesh(dim, nv, ne, nb, sdim);

    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            mesh.AddVertex(i%3, j%3);
        }
    }

    mesh.AddQuad(0,1,4,3);
    mesh.AddQuad(1,2,5,4);
    mesh.AddQuad(3,4,7,6);
    mesh.AddQuad(4,5,8,7);
    mesh.FinalizeQuadMesh();

    mesh.Save("mesh.mesh");
}