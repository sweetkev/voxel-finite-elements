#include "mfem.hpp"
#include <tuple>
#include <string>

using namespace mfem;

std::tuple<std::unordered_map<std::string, int>, Mesh> coarsenMesh(int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex);
std::unordered_map<std::string, int> addCoarseVertices(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex);
bool pixelNearby(int i, int j, std::unordered_map<std::string, int> coord_to_fine_vertex);
void addQuads(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex, std::unordered_map<std::string, int> coord_to_coarse_vertex);