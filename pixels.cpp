#include "make_mesh.hpp"
#include "ppm.hpp"
#include "../mfem/mfem-4.8/mfem.hpp"
#include <fstream>

int main(int argc, char *argv[]) {
    //parse options
    std::string pgm_file = "pgm_files/small.pgm";

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&pgm_file, "-f", "--file", "pgm file to use");
    args.ParseCheck();

    makeMesh(pgm_file);
}