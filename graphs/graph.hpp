#pragma once

#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>

using namespace mfem;

class Graph
{
public:
    /**
     * Returns the graph where nodes represent cells, and edges exist
     *  between cells that share DoFs
     */
    Graph(const FiniteElementSpace &fes, const PixelImage &image_);

    /** Returns the mfem table representing the graph */
    const Table &GetGraph() { return graph; }

    /** Returns index of the connection between node i and node j.
     * If there is no connection between node i and node j established in the
     * table, then the return value is -1. */
    int operator() (int i, int j) const { return graph(i,j); }

    /** Returns the number of nodes in the graph */
    int Size() const { return graph.Size(); }

    /** Returns a graph over a coarsened mesh */
    Graph CoarsenGraph();

private:
    Graph(Table &graph_,
          std::unordered_map<Coord, std::vector<int>> &coord_to_node_,
          std::vector<Coord> &node_to_coord_);

    Table graph;

    // Maps between coordiates and cells
    std::unordered_map<Coord, std::vector<int>> coord_to_node;
    std::vector<Coord> node_to_coord;
    std::vector<Coord> grid_cells; // Occupied grid cells

    /**
     * Labels each node with the index of the node it
     *  corresponds to in the coarse graph. It also creates the
     *  cell-to-node map for the coarse graph.
    */
    Array<int> LabelGraph(std::unordered_map<Coord, std::vector<int>> &coarse_coord_to_node,
                          std::vector<Coord> &coarse_node_to_coord);

    /** Builds coarse graph according to graph labeling */
    Table BuildCoarseGraph(Array<int> &graph_labeling,
        int coarse_ne);
};
