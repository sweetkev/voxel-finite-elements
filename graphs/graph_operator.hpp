#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include "graph.hpp"
#include <unordered_map>

using namespace mfem;



class GraphOperator
{
public:

private:
    Mesh reference_mesh;
    SparseMatrix Aref;
    SparseMatrix V;

    void CreateReferenceMesh(int dim);
};