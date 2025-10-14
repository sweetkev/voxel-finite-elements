#include "mfem.hpp"
#include "coarsen_mesh.hpp"
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <string>

using namespace mfem;

std::tuple<std::unordered_map<std::string, int>, Mesh> coarsenMesh(Mesh &fine_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex) {

    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //map that takes a coordinate pair to a vertex number
    std::unordered_map<std::string, int> coord_to_coarse_vertex = findCoarseVertices(fine_mesh, m, n, coord_to_fine_vertex);

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

std::unordered_map<std::string, int> findCoarseVertices(Mesh &fine_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex) {
    int p = std::ceil(0.5*m), q = std::ceil(0.5*n);

    std::unordered_map<std::string, int> coord_to_coarse_vertex;

    for(int j = 0; j < q; j++) {
        for(int i = 0; i < p; i++) {
            std::string coord = std::to_string(2*i) + " " + std::to_string(2*j);
        }
    }

}

int numCoarseElements(Mesh &fine_mesh, int m, int n) {

}