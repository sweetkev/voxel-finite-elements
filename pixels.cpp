#include "ppm.hpp"
#include "mfem.hpp"
#include "pixel_mesh.hpp"
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

    //make mesh from file
    PixelMesh fine_mesh(pgm_file);

    //coarsen mesh
    PixelMesh coarse_mesh = fine_mesh.CoarsenMesh();
}