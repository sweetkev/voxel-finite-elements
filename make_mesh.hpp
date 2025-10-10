#ifndef MAKE_MESH_HPP
#define MAKE_MESH_HPP

#include "../mfem/mfem-4.8/mfem.hpp"
#include "ppm.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;

bool adjacentPixelFilled(int i, int j, PixelImage image);
int numElements(PixelImage image);
void makeMesh(std::string pgm_file);

#endif