#ifndef MAKE_MESH_HPP
#define MAKE_MESH_HPP

#include "mfem.hpp"
#include "ppm.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <tuple>

using namespace mfem;

std::tuple<std::unordered_map<std::string, int>, Mesh> makeMesh(PixelImage image);
bool adjacentPixelFilled(int i, int j, PixelImage image);
int numElements(PixelImage image);
std::unordered_map<std::string, int> findVertices(PixelImage image);

#endif