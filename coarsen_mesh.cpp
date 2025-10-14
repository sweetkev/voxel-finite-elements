#include "mfem.hpp"
#include "coarsen_mesh.hpp"
#include <iostream>
#include <tuple>
#include <string>

using namespace mfem;

std::tuple<std::unordered_map<std::string, int>, Mesh> coarsenMesh(Mesh *fine_mesh, int m, int n) {

    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //map that takes a coordinate pair to a vertex number
    std::unordered_map<std::string, int> coord_to_coarse_vertex = findCoarseVertices(fine_mesh, m, n);

    //number of vertices
    int nv = coord_to_coarse_vertex.size();

    //find number of elements
    int ne = numCoarseElements(fine_mesh, m, n);

    //ignoring boundary for now
    int nb = 0;

    //initialize mesh
    Mesh mesh(dim, nv, ne, nb, sdim);


    Mesh coarse_mesh;

    return std::make_tuple(coord_to_coarse_vertex, coarse_mesh);
}

std::unordered_map<std::string, int> findCoarseVertices(Mesh *fine_mesh, int m, int n) {

}

int numCoarseElements(Mesh *fine_mesh, int m, int n) {

}