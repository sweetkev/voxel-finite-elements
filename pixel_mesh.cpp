#include "ppm.hpp"
#include "mfem.hpp"
#include "pixel_mesh.hpp"
#include <iostream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

PixelMesh::PixelMesh(std::string pgm_file) {
    PixelImage image(pgm_file);
    std::tie(coord_to_vertex, mesh) = MakeMesh(image);
    width = image.Width();
    height = image.Height();
}

PixelMesh::PixelMesh(Mesh &mesh_, std::unordered_map<std::string, int> coord_to_vertex_, int width_, int height_) 
    : mesh(mesh_), coord_to_vertex(coord_to_vertex_), width(width_), height(height_) 
{}

PixelMesh PixelMesh::CoarsenMesh() {
    int m = width, n = height;
    std::unordered_map<std::string, int> coord_to_fine_vertex = coord_to_vertex;

    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //number of vertices, elements, and boundary elements. Will be allocated by FinalizeMesh() later
    int nv = 0;
    int ne = 0;
    int nb = 0;

    //initialize mesh
    Mesh coarse_mesh(dim, nv, ne, nb, sdim);

    //add vertices
    std::unordered_map<std::string, int> coord_to_coarse_vertex = AddCoarseVertices(coarse_mesh, m, n, coord_to_fine_vertex);

    //add quads
    AddCoarseQuads(coarse_mesh, m, n, coord_to_fine_vertex, coord_to_coarse_vertex);

    coarse_mesh.FinalizeMesh();
    coarse_mesh.Save("coarse_mesh.mesh");

    PixelMesh c_mesh(coarse_mesh, coord_to_coarse_vertex, m, n);
    return c_mesh;
}

SparseMatrix PixelMesh::CreateProlongation(PixelMesh fine_mesh) {
    std::unordered_map<std::string, int> fine_coord_to_vertex = fine_mesh.GetVertexMap();
    
    //initialize matrix
    int nrows = fine_coord_to_vertex.size();
    int ncols = coord_to_vertex.size();
    SparseMatrix P(nrows,ncols);

    //add entries
    int p = std::ceil(0.5*(width))+1, q = std::ceil(0.5*(height))+1;

    int vertex = 0;
    for(int j = 0; j < q; j++) {
        for(int i = 0; i < p; i++) { 
            int x = 2*i, y = 2*(q-j-1);
            std::string coord = std::to_string(x) + " " + std::to_string(y);
            if(coord_to_vertex.count(coord) != 0) {
                AddVertexInfluence(x,y,P,fine_coord_to_vertex);
            }
        }
    }

    //finalize and return
    P.Finalize();
    return P;
}

std::tuple<std::unordered_map<std::string, int>, Mesh> PixelMesh::MakeMesh(PixelImage image) { 
    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //number of vertices, elements, and boundary elements. Will be allocated by FinalizeMesh() later
    int nv = 0;
    int ne = 0;
    int nb = 0;

    //initialize mesh
    Mesh fine_mesh(dim, nv, ne, nb, sdim);

    //add vertices
    int m = image.Width(), n = image.Height();
    std::unordered_map<std::string, int> coord_to_vertex = AddVertices(image, fine_mesh, m, n);

    //add quads
    AddQuads(image, fine_mesh, m, n, coord_to_vertex);

    fine_mesh.FinalizeMesh();
    fine_mesh.Save("fine_mesh.mesh");

    return std::make_tuple(coord_to_vertex, fine_mesh);
}

bool PixelMesh::AdjacentPixelFilled(int i, int j, PixelImage image) {
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

std::unordered_map<std::string, int> PixelMesh::AddVertices(PixelImage image, Mesh &mesh, int m, int n) {
    std::unordered_map<std::string, int> coord_to_vertex;

    int vertex = 0;
    for(int j = 0; j < n+1; j++) {
        for(int i = 0; i < m+1; i++) {
            if(AdjacentPixelFilled(i,j,image)) {
                int x = i, y = n-j;
                mesh.AddVertex(x,y);
                std::string coord = std::to_string(x) + " " + std::to_string(y);
                coord_to_vertex.insert({coord, vertex});
                vertex++;
            }
        }
    }

    return coord_to_vertex;
}

void PixelMesh::AddQuads(PixelImage image, Mesh &mesh, int m, int n, std::unordered_map<std::string, int> coord_to_vertex) {
    for(int j = 0; j < n; j++) {
        for(int i = 0; i < m; i++) {
            //add quad for filled in pixels
            if(image.operator()(i,j) != 0) {
                int x = i, y = n-j;
                int v1 = coord_to_vertex[std::to_string(x) + " " + std::to_string(y)];
                int v2 = coord_to_vertex[std::to_string(x) + " " + std::to_string(y-1)];
                int v3 = coord_to_vertex[std::to_string(x+1) + " " + std::to_string(y-1)];
                int v4 = coord_to_vertex[std::to_string(x+1) + " " + std::to_string(y)];
                mesh.AddQuad(v1,v2,v3,v4);
            }
        }
    }
}

bool PixelMesh::PixelNearby(int x, int y, std::unordered_map<std::string, int> coord_to_fine_vertex) {
    std::string coord1 = std::to_string(x - 1) + " " + std::to_string(y - 1);
    std::string coord2 = std::to_string(x + 1) + " " + std::to_string(y - 1);
    std::string coord3 = std::to_string(x + 1) + " " + std::to_string(y + 1);
    std::string coord4 = std::to_string(x - 1) + " " + std::to_string(y + 1);
    if((coord_to_fine_vertex.count(coord1) != 0) || (coord_to_fine_vertex.count(coord2) != 0) || (coord_to_fine_vertex.count(coord3) != 0) || (coord_to_fine_vertex.count(coord4) != 0)) {
        return true;
    }

    return false;
}

std::unordered_map<std::string, int> PixelMesh::AddCoarseVertices(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex) {
    int p = std::ceil(0.5*(m))+1, q = std::ceil(0.5*(n))+1;

    std::unordered_map<std::string, int> coord_to_coarse_vertex;

    int vertex = 0;
    for(int j = 0; j < q; j++) {
        for(int i = 0; i < p; i++) {
            int x = 2*i, y = 2*(q-j-1);
            std::string coord = std::to_string(x) + " " + std::to_string(y);
            if(PixelNearby(x,y,coord_to_fine_vertex)) {
                coarse_mesh.AddVertex(x,y);
                coord_to_coarse_vertex.insert({coord,vertex});
                vertex++;
            }
        }
    }

    return coord_to_coarse_vertex;
}

void PixelMesh::AddCoarseQuads(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex, std::unordered_map<std::string, int> coord_to_coarse_vertex) {
    int p = std::ceil(0.5*m), q = std::ceil(0.5*n);
    
    for(int j = 0; j < q; j++) {
        for(int i = 0; i < p; i++) {
            //add quad for filled in pixels
            int x = 2*i, y = 2*(q-j-1);
            std::string coord = std::to_string(x+1) + " " + std::to_string(y-1);
            if(coord_to_fine_vertex.count(coord) != 0) {
                int v1 = coord_to_coarse_vertex[std::to_string(x) + " " + std::to_string(y)];
                int v2 = coord_to_coarse_vertex[std::to_string(x) + " " + std::to_string(y-2)];
                int v3 = coord_to_coarse_vertex[std::to_string(x+2) + " " + std::to_string(y-2)];
                int v4 = coord_to_coarse_vertex[std::to_string(x+2) + " " + std::to_string(y)];
                coarse_mesh.AddQuad(v1,v2,v3,v4);
            }
        }
    }
}

void PixelMesh::AddVertexInfluence(int x, int y, SparseMatrix &P, std::unordered_map<std::string, int> fine_coord_to_vertex) {
    std::string coarse_coord = std::to_string(x) + " " + std::to_string(y);
    int i,j;
    j = coord_to_vertex[coarse_coord];
    
    std::string coord;

    //current (middle middle) vertex
    coord = coarse_coord;
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,1);
    }
    
    //bottom left vertex
    coord = std::to_string(x-1) + " " + std::to_string(y-1);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.25);
    }

    //bottom middle vertex
    coord = std::to_string(x) + " " + std::to_string(y-1);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.5);
    }

    //bottom right vertex
    coord = std::to_string(x+1) + " " + std::to_string(y-1);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.25);
    }

    //middle left vertex
    coord = std::to_string(x-1) + " " + std::to_string(y);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.5);
    }

    //middle right vertex
    coord = std::to_string(x+1) + " " + std::to_string(y);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.5);
    }

    //top left vertex
    coord = std::to_string(x-1) + " " + std::to_string(y+1);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.25);
    }

    //top middle vertex
    coord = std::to_string(x) + " " + std::to_string(y+1);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.5);
    }

    //top right vertex
    coord = std::to_string(x+1) + " " + std::to_string(y+1);
    if(fine_coord_to_vertex.count(coord) != 0) {
        i = fine_coord_to_vertex[coord];
        P.Add(i,j,0.25);
    }
}