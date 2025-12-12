#include "sparsify/DSB.hpp"
#include "graph/Graph.hpp"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace sparsify
{
    struct UEdge
    {
        int    u;
        int    v;
        double w;
    };

    /* keep edges (u, v) for which v > u: each undirected edge is sampled once */
    static std::vector<UEdge> collect_unique_edges (const Graph& G)
    {
        
    }
}