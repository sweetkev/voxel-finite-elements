#include "make_mesh.hpp"
#include "ppm.hpp"
//ensure path to mfem.hpp is correct
#include "../mfem/mfem-4.8/mfem.hpp"
#include <iostream>
#include <string>

using namespace mfem;

void makeMesh(std::string pgm_file) {
    //make PixelImage
    PixelImage image(pgm_file);
    int m = image.Width(), n = image.Height(); 
    
    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //find number of vertices, and start of index in each row stored in an array
    //the entry in r_start[j] gives the index of the first vertex in row j
    int r_start[n+1];                                                   
    int nv = 0;
    for(int j = 0; j < n+1; j++) {
        r_start[j] = nv;
        for(int i = 0; i < m+1; i++) {
            if(adjacentPixelFilled(i,j,image)) {
                nv++;
            }
        }
    }
    
    //find number of elements
    int ne = numElements(image);

    //ignoring boundary for now
    int nb = 0;

    //initialize mesh
    Mesh mesh(dim, nv, ne, nb, sdim);

    //add vertices (uses integer coordinates, so each pixel will be 1x1)
    for(int j = 0; j < n+1; j++) {
        for(int i = 0; i < m+1; i++) {
            if(adjacentPixelFilled(i,j,image)) {
                mesh.AddVertex(i,n-j);
            }
        }
    }

    //add quads
    int ui = 0;
    int di = 0;
    for(int j = 0; j < n; j++) {
        //"upper index" and "down index". Keeps track of vertex index in rows of vertices above and below pixels.
        ui = 0;
        di = 0;
        for(int i = 0; i < m; i++) {
            //add quad for filled in pixels
            if(image.operator()(i,j) != 0) {
                mesh.AddQuad(r_start[j+1]+di, r_start[j+1]+di+1, r_start[j]+ui+1, r_start[j]+ui);
                ui++;
                di++;
                continue;
            }

            //if current pixel is not filled in, we progress down and upper indices if filled in pixels are above or below or left
            if(i != 0) {
                if(image.operator()(i-1,j) != 0) {
                    ui++;
                    di++;
                    continue;
                }
            }
            if(j != 0) {
                if(image.operator()(i,j-1) != 0) {
                    ui++;
                }
            }
            if(j != n-1) {
                if(image.operator()(i,j+1) != 0) {
                    di++;
                }
            }
        }
    }

    mesh.Save("mesh.mesh");
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

