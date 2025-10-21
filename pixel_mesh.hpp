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

        PixelMesh CoarsenMesh();

        Mesh& GetMesh() { return mesh; }
        std::unordered_map<std::string, int> GetVertexMap() { return coord_to_vertex; }


    private:
        Mesh mesh;
        std::unordered_map<std::string, int> coord_to_vertex;
        int width = -1;
        int height = -1;

        std::tuple<std::unordered_map<std::string, int>, Mesh> MakeMesh(PixelImage image);

        bool AdjacentPixelFilled(int i, int j, PixelImage image);
        std::unordered_map<std::string, int> AddVertices(PixelImage image, Mesh &mesh, int m, int n);
        void AddQuads(PixelImage image, Mesh &mesh, int m, int n, std::unordered_map<std::string, int> coord_to_vertex);

        bool PixelNearby(int i, int j, std::unordered_map<std::string, int> coord_to_fine_vertex);
        std::unordered_map<std::string, int> AddCoarseVertices(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex);
        void AddCoarseQuads(Mesh &coarse_mesh, int m, int n, std::unordered_map<std::string, int> coord_to_fine_vertex, std::unordered_map<std::string, int> coord_to_coarse_vertex);
};