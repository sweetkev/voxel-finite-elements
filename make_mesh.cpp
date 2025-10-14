#include "make_mesh.hpp"
#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

std::tuple<std::unordered_map<std::string, int>, Mesh> makeMesh(PixelImage image) {
    int m = image.Width(), n = image.Height(); 
    
    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //map that takes a coordinate pair to a vertex number
    std::unordered_map<std::string, int> coord_to_vertex = findVertices(image);

    //number of vertices
    int nv = coord_to_vertex.size();

    //find number of elements
    int ne = numElements(image);

    //ignoring boundary for now
    int nb = 0;

    //initialize mesh
    Mesh mesh(dim, nv, ne, nb, sdim);

    //add vertices
    addVertices(mesh, m, n, coord_to_vertex);

    //add quads
    addQuads(image, mesh, m, n, coord_to_vertex);

    mesh.FinalizeMesh();
    mesh.Save("fine_mesh.mesh");

    return std::make_tuple(coord_to_vertex, mesh);
}

/*
returns true when a pixel adjacent to vertex is filled in
*/
bool adjacentPixelFilled(int i, int j, PixelImage image) {
    int m = image.Width(), n = image.Height();

    //bottom row cases
    if(j == n) {
        //bottom left vertex
        if(i == 0) {
            return image.operator()(i,j-1) != 0;
        }
        //bottom right vertex
        else if(i == m) {
            return image.operator()(i-1,j-1) != 0;
        }
        //other bottom vertices
        else {
            return (image.operator()(i-1,j-1) != 0) || (image.operator()(i,j-1) != 0);
        }
    }

    //top row cases
    if(j == 0) {
        //top left vertex
        if(i == 0) {
            return image.operator()(i,j) != 0;
        }
        //top right vertex
        else if(i == m) {
            return image.operator()(i-1,j) != 0;
        }
        //other top vertices
        else {
            return (image.operator()(i,j) != 0) || (image.operator()(i-1,j) != 0);
        }
    }

    //all other rows
    //left vertices
    if(i == 0) {
        return (image.operator()(i,j-1) != 0) || (image.operator()(i,j) != 0);
    }
    //right vertices
    else if(i == m) {
        return (image.operator()(i-1,j-1) != 0) || (image.operator()(i-1,j) != 0);
    }
    //all other vertices (interior vertices)
    else {
        return (image.operator()(i-1,j-1) != 0) || (image.operator()(i-1,j) != 0) || (image.operator()(i,j-1) != 0) || (image.operator()(i,j) != 0);
    }
}

/*
returns number of elements in mesh generated from image
*/
int numElements(PixelImage image) {
    int m = image.Width(), n = image.Height();

    int ne = 0;
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            if(image.operator()(i,j) != 0) {
                ne++;
            }
        }
    }

    return ne;
}

/*
returns an unordered map that maps coordinate values to vertices
*/
std::unordered_map<std::string, int> findVertices(PixelImage image) {
    int m = image.Width(), n = image.Height();
    std::unordered_map<std::string, int> coord_to_vertex;
    int vertex = 0;

    for(int j = 0; j < n+1; j++) {
        for(int i = 0; i < m+1; i++) {
            if(adjacentPixelFilled(i,j,image)) {
                std::string coord = std::to_string(i) + " " + std::to_string(j);
                coord_to_vertex.insert({coord, vertex});
                vertex++;
            }
        }
    }

    return coord_to_vertex;
}

/*
adds vertices to mesh at integer indicies, so pixels will be 1x1
*/
void addVertices(Mesh &mesh, int m, int n, std::unordered_map<std::string, int> coord_to_vertex) {
    for(int j = 0; j < n+1; j++) {
        for(int i = 0; i < m+1; i++) {
            std::string coord = std::to_string(i) + " " + std::to_string(j);
            if(coord_to_vertex.count(coord) != 0) {
                mesh.AddVertex(i,n-j);
            }
        }
    }
}

/*
adds quads for filled pixels
*/
void addQuads(PixelImage image, Mesh &mesh, int m, int n, std::unordered_map<std::string, int> coord_to_vertex) {
    for(int j = 0; j < n; j++) {
        for(int i = 0; i < m; i++) {
            //add quad for filled in pixels
            if(image.operator()(i,j) != 0) {
                int v1 = coord_to_vertex[std::to_string(i) + " " + std::to_string(j)];
                int v2 = coord_to_vertex[std::to_string(i) + " " + std::to_string(j+1)];
                int v3 = coord_to_vertex[std::to_string(i+1) + " " + std::to_string(j+1)];
                int v4 = coord_to_vertex[std::to_string(i+1) + " " + std::to_string(j)];
                mesh.AddQuad(v1,v2,v3,v4);
            }
        }
    }
}