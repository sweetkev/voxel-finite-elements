#include "mfem.hpp"
#include <tuple>
#include <string>

std::tuple<std::unordered_map<std::string, int>, Mesh> coarsenMesh(Mesh *fine_mesh, int m, int n);
std::unordered_map<std::string, int> findCoarseVertices(Mesh *fine_mesh, int m, int n);
int numCoarseElements(Mesh *fine_mesh, int m, int n);