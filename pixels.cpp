#include "make_mesh.hpp"
#include "coarsen_mesh.hpp"
#include "ppm.hpp"
#include "mfem.hpp"
#include <fstream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

int main(int argc, char *argv[]) {
    //parse options
    std::string pgm_file = "pgm_files/small.pgm";

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&pgm_file, "-f", "--file", "pgm file to use");
    args.ParseCheck();

    //make PixelImage
    PixelImage image(pgm_file);

    //make fine mesh from image
    std::unordered_map<std::string, int> coord_to_fine_vertex;
    Mesh fine_mesh;
    std::tie(coord_to_fine_vertex, fine_mesh) = makeMesh(image);

    //coarsen mesh
    std::unordered_map<std::string, int> coord_to_coarse_vertex;
    Mesh coarse_mesh;
    std::tie(coord_to_coarse_vertex, coarse_mesh) = coarsenMesh(image.Width(), image.Height(), coord_to_fine_vertex);
}