#include "mfem.hpp"

using namespace mfem;

class Graph 
{
    public:
        /** 
         * Returns the graph where nodes represent elements, and edges exist
         *  between elements that share DoFs
         */
        Graph(const FiniteElementSpace &fes);

        const SparseMatrix &GetGraph() { return graph; }

    private:
        SparseMatrix graph;
        Array<int> node_to_element;

        /** Returns true if the two arrays of DoFs share any DoFs */
        bool SharesDof(const Array<int> &idofs, const Array<int> &jdofs);
};