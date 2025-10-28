#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

class PixelMesh
{
    public:
        PixelMesh(std::string pgm_file);
        PixelMesh(Mesh &mesh_, std::unordered_map<std::string, int> coord_to_vertex_, int width_, int height_);

        /**
            Creates coarse mesh from current mesh
        */
        PixelMesh CoarsenMesh();

        Mesh& GetMesh() { return mesh; }
        std::unordered_map<std::string, int> GetVertexMap() { return coord_to_vertex; }


    private:
        Mesh mesh;
        std::unordered_map<std::string, int> coord_to_vertex;
        int width = -1;
        int height = -1;

        /**
            Makes mesh given PixelImage, with elements given by filled-in pixels
        */
        std::tuple<std::unordered_map<std::string, int>, Mesh> MakeMesh(PixelImage image);

        /**
            Returns true when a pixel adjacent to vertex is filled in. Else returns false.
        */
        bool AdjacentPixelFilled(int i, int j, PixelImage image);

        /**
            Adds vertices to mesh at integer indicies, so pixels will be 1x1. Returns map taking coordinates to vertex number.
        */
        std::unordered_map<std::string, int> AddVertices(PixelImage image, Mesh &mesh, int m, int n);
        
        /**
            Adds quads for filled pixels
        */
        void AddQuads(PixelImage image, Mesh &mesh, int m, int n, std::unordered_map<std::string, int> coord_to_vertex);

        /**
            Returns true if any adjacent coarse pixel should be filled in. Else returns false
        */
        bool PixelNearby(int i, int j, std::unordered_map<std::string, int> coord_to_fine_vertex);
        
        /**
            Returns unorderd map that takes coordinates to coarse vertex, and adds vertices to coarse mesh
        */
        std::unordered_map<std::string, int> AddCoarseVertices(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex);
        
        /**
            Adds quads to coarse mesh
        */
        void AddCoarseQuads(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex, std::unordered_map<std::string, int> coord_to_coarse_vertex);
};