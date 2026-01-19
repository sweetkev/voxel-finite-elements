#include "mfem.hpp"
#include "ppm.hpp"
#include "pixel_mesh.hpp"
#include <unordered_map>

using namespace mfem;

class Graph 
{
public:
    Graph(Table &graph_,
        Array<int> &node_to_element_,
        std::vector<std::vector<int>> &element_to_node_,
        std::unordered_map<Coord, int> &coord_to_element_,
        std::vector<Coord> &element_to_coord_,
        const PixelImage &image_)
        : graph(graph_),
        node_to_element(node_to_element_),
        element_to_node(element_to_node_),
        coord_to_element(coord_to_element_),
        element_to_coord(element_to_coord_),
        image(image_) { }

    /** 
     * Returns the graph where nodes represent elements, and edges exist
     *  between elements that share DoFs
     */
    Graph(const FiniteElementSpace &fes, const PixelImage &image_);

    /** Returns the mfem table representing the graph */
    const Table &GetGraph() { return graph; }

    /** Returns the geometric mesh element represented by node i */
    int GetNodeElement(int i) { return node_to_element[i]; }

    //** Returns the array of nodes that 'occupy' element i */
    std::vector<int> GetElementNodes(int i) 
        { return element_to_node[i]; }

    /** Returns the bottom-left coordinates of element i */
    Coord GetElementCoord(int i) const { return element_to_coord[i]; }

    /** 
     * Returns the element index of the element with bottom-left 
     *  coordinates coord 
     */
    int GetElementIndex(Coord coord) const
    {
        auto it = coord_to_element.find(coord);
        if (it != coord_to_element.end())
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

    // Maps between nodes and geometric elements
    Array<int> node_to_element;
    std::vector<std::vector<int>> element_to_node;

    // Maps between coordiates and elements
    std::unordered_map<Coord, int> coord_to_element;
    std::vector<Coord> element_to_coord;

    /** Returns true if the two arrays of DoFs share any DoFs */
    bool SharesDof(const Array<int> &idofs, const Array<int> &jdofs);

    /** Creates maps between coordinates and elements */
    void CreateCoordElementMaps(
        std::unordered_map<Coord, int> &coord_to_vertex,
        std::vector<Coord> &vertex_to_coord, const PixelImage &image_);

    /** 
     * Labels each node with the index of the node it 
     *  corresponds to in the coarse graph. It also creates the 
     *  element-to-node map for the coarse graph.
    */
    Array<int> LabelGraph(
        Array<int> &coarse_node_to_element, 
        std::vector<std::vector<int>> &coarse_element_to_node,
        const std::unordered_map<Coord, int> &coarse_coord_to_element,
        const std::vector<Coord> &coarse_element_to_coord, 
        const PixelImage &coarse_image);

    /** Builds coarse graph according to graph labeling */
    Table BuildCoarseGraph(Array<int> &graph_labeling,
        int coarse_ne);
    
    /** Removes diagonal entries from graph */
    void RemoveSelfConnections(Table &coarse_graph);
};