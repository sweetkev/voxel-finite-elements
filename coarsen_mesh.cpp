#include "mfem.hpp"
#include "coarsen_mesh.hpp"
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <string>

using namespace mfem;

std::tuple<std::unordered_map<std::string, int>, Mesh> coarsenMesh(int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex) {

    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //number of vertices, elements, and boundary elements. Will be allocated by FinalizeMesh() later
    int nv = 0;
    int ne = 0;
    int nb = 0;

    //initialize mesh
    Mesh coarse_mesh(dim, nv, ne, nb, sdim);

    //add vertices
    std::unordered_map<std::string, int> coord_to_coarse_vertex = addCoarseVertices(coarse_mesh, m, n, coord_to_fine_vertex);

    //add quads
    addQuads(coarse_mesh, m, n, coord_to_fine_vertex, coord_to_coarse_vertex);

    coarse_mesh.FinalizeMesh();
    coarse_mesh.Save("coarse_mesh.mesh");

    return std::make_tuple(coord_to_coarse_vertex, coarse_mesh);
}

/*
returns unorderd map that takes coordinates to coarse vertex, and adds vertices to coarse mesh
*/
std::unordered_map<std::string, int> addCoarseVertices(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex) {
    int p = std::ceil(0.5*(m+1)), q = std::ceil(0.5*(n+1));

    std::unordered_map<std::string, int> coord_to_coarse_vertex;

    int vertex = 0;
    for(int j = 0; j < q; j++) {
        for(int i = 0; i < p; i++) {
            std::string coord = std::to_string(2*i) + " " + std::to_string(2*j);
            if(pixelNearby(i,j,coord_to_fine_vertex)) {
                coarse_mesh.AddVertex(2*i,(2*(q-j-1)));
                coord_to_coarse_vertex.insert({coord,vertex});
                vertex++;
            }
        }
    }

    return coord_to_coarse_vertex;

}

/*
returns true if any adjacent coarse pixel should be filled in. Else returns false
*/
bool pixelNearby(int i, int j, std::unordered_map<std::string, int> coord_to_fine_vertex) {
    std::string coord1 = std::to_string(2*i - 1) + " " + std::to_string(2*j - 1);
    std::string coord2 = std::to_string(2*i + 1) + " " + std::to_string(2*j - 1);
    std::string coord3 = std::to_string(2*i + 1) + " " + std::to_string(2*j + 1);
    std::string coord4 = std::to_string(2*i - 1) + " " + std::to_string(2*j + 1);
    if((coord_to_fine_vertex.count(coord1) != 0) || (coord_to_fine_vertex.count(coord2) != 0) || (coord_to_fine_vertex.count(coord3) != 0) || (coord_to_fine_vertex.count(coord4) != 0)) {
        return true;
    }

    return false;
}

/*
adds quads to coarse mesh
*/
void addQuads(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex, std::unordered_map<std::string, int> coord_to_coarse_vertex) {
    int p = std::ceil(0.5*m), q = std::ceil(0.5*n);
    
    for(int j = 0; j < q; j++) {
        for(int i = 0; i < p; i++) {
            //add quad for filled in pixels
            std::string coord = std::to_string(2*i+1) + " " + std::to_string(2*j+1);
            if(coord_to_fine_vertex.count(coord) != 0) {
                int v1 = coord_to_coarse_vertex[std::to_string(2*i) + " " + std::to_string(2*j)];
                int v2 = coord_to_coarse_vertex[std::to_string(2*i) + " " + std::to_string(2*(j+1))];
                int v3 = coord_to_coarse_vertex[std::to_string(2*(i+1)) + " " + std::to_string(2*(j+1))];
                int v4 = coord_to_coarse_vertex[std::to_string(2*(i+1)) + " " + std::to_string(2*j)];
                coarse_mesh.AddQuad(v1,v2,v3,v4);
            }
        }
    }
}