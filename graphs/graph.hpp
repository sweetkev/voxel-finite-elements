#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>

using namespace mfem;

class Graph 
{
public:
    Graph(Table &graph_,
        Array<int> &node_to_cell_,
        std::vector<std::vector<int>> &cell_to_node_,
        std::unordered_map<Coord, int> &coord_to_cell_,
        std::vector<Coord> &cell_to_coord_,
        const PixelImage &image_)
        : graph(graph_),
        node_to_cell(node_to_cell_),
        cell_to_node(cell_to_node_),
        coord_to_cell(coord_to_cell_),
        cell_to_coord(cell_to_coord_),
        image(image_) { }

    /** 
     * Returns the graph where nodes represent cells, and edges exist
     *  between cells that share DoFs
     */
    Graph(const FiniteElementSpace &fes, const PixelImage &image_);

    /** Returns the mfem table representing the graph */
    const Table &GetGraph() { return graph; }

    /** Returns the geometric mesh cell represented by node i */
    int GetNodeCell(int i) { return node_to_cell[i]; }

    //** Returns the array of nodes that 'occupy' cell i */
    std::vector<int> GetCellNodes(int i) 
        { return cell_to_node[i]; }

    /** Returns the bottom-left coordinates of cell i */
    Coord GetCellCoord(int i) const { return cell_to_coord[i]; }

    /** 
     * Returns the cell index of the cell with bottom-left 
     *  coordinates coord 
     */
    int GetCellIndex(Coord coord) const
    {
        auto it = coord_to_cell.find(coord);
        if (it != coord_to_cell.end())
        {
            return it->second;
        }
        else
        {
            return -1;
        }
    }

    /** Returns a graph over a coarsened mesh */
    Graph CoarsenGraph();

private:
    Table graph;

    PixelImage image;

    // Maps between nodes and geometric cells
    Array<int> node_to_cell;
    std::vector<std::vector<int>> cell_to_node;

    // Maps between coordiates and cells
    std::unordered_map<Coord, int> coord_to_cell;
    std::vector<Coord> cell_to_coord;

    /** Creates maps between coordinates and cells */
    void CreateCoordCellMaps(
        std::unordered_map<Coord, int> &coord_to_vertex,
        std::vector<Coord> &vertex_to_coord, const PixelImage &image_);

    /** 
     * Labels each node with the index of the node it 
     *  corresponds to in the coarse graph. It also creates the 
     *  cell-to-node map for the coarse graph.
    */
    Array<int> LabelGraph(
        Array<int> &coarse_node_to_cell, 
        std::vector<std::vector<int>> &coarse_cell_to_node,
        const std::unordered_map<Coord, int> &coarse_coord_to_cell,
        const std::vector<Coord> &coarse_cell_to_coord, 
        const PixelImage &coarse_image);

    /** Builds coarse graph according to graph labeling */
    Table BuildCoarseGraph(Array<int> &graph_labeling,
        int coarse_ne);
};