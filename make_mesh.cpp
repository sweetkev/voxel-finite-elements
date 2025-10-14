#include "make_mesh.hpp"
#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace mfem;

std::unordered_map<std::string, int> makeMesh(std::string pgm_file) {
    //make PixelImage
    PixelImage image(pgm_file);
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

    //add vertices (uses integer coordinates, so each pixel will be 1x1)
    for(int j = 0; j < n+1; j++) {
        for(int i = 0; i < m+1; i++) {
            std::string coord = std::to_string(i) + " " + std::to_string(j);
            if(coord_to_vertex.count(coord) != 0) {
                mesh.AddVertex(i,n-j);
            }
        }
    }

    //add quads
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

    mesh.Finalize();
    mesh.Save("fine_mesh.mesh");

    return coord_to_vertex;
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