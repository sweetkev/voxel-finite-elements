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

    /** Returns the array of nodes connected to current node */
    const Array<int> GetConnectedNodes(int node) const;

    /** Returns the elements with coordinate Coord */
    const std::vector<int> &GetCoordElements(Coord coord) const { return coord_to_node.at(coord); }

    /** Returns the coord of element e */
    const Coord &GetElementCoord(int e) const { return node_to_coord[e]; }

    /** Returns the occupied grid cells */
    const std::vector<Coord> &GetGridCells() const { return grid_cells; }

    /** Returns index of the connection between node i and node j.
     * If there is no connection between node i and node j established in the
     * table, then the return value is -1. */
    int operator() (int i, int j) const { return graph(i,j); }

    /** Returns the number of nodes in the graph */
    int Size() const { return graph.Size(); }

    /** Returns a graph over a coarsened mesh */
    Graph CoarsenGraph();

    /** Returns the fine-to-coarse element labeling for the coarse graph */
    const Array<int> &GetGraphLabeling() const { return graph_labeling; }

private:
    Graph(Table &graph_,
          std::unordered_map<Coord, std::vector<int>> &coord_to_node_,
          std::vector<Coord> &node_to_coord_, std::vector<Coord> &grid_cells_,
          const Array<int> &graph_labeling_ = Array<int>());

    Table graph;

    // Maps between coordiates and nodes
    std::unordered_map<Coord, std::vector<int>> coord_to_node;
    std::vector<Coord> node_to_coord;
    std::vector<Coord> grid_cells; // Occupied grid cells
    Array<int> graph_labeling;

    /**
     * Labels each node with the index of the node it
     *  corresponds to in the coarse graph. It also creates the
     *  cell-to-node map for the coarse graph.
    */
    Array<int> LabelGraph(std::unordered_map<Coord, std::vector<int>> &coarse_coord_to_node,
                          std::vector<Coord> &coarse_node_to_coord,
                          std::vector<Coord> &coarse_grid_cells);

    /** Builds coarse graph according to graph labeling */
    Table BuildCoarseGraph(const Array<int> &graph_labeling,
        int coarse_ne);
};
