#include "graph/Graph.hpp"
#include <iostream>

Graph::Graph(int n,
             std::vector<Edge>&& edges,
             std::vector<int>&& row_ptr,
             std::vector<int>&& col_idx,
             std::vector<double>&& adj_w,
             std::vector<double>&& degree,
             bool connected)
  : n_(n),
    m_(static_cast<int>(edges.size())),
    edges_(std::move(edges)),
    row_ptr_(std::move(row_ptr)),
    col_idx_(std::move(col_idx)),
    adj_w_(std::move(adj_w)),
    degree_(std::move(degree)),
    connected_(connected) {}



void Graph::build_from_edges ()
{   
    if (built_)         throw std::runtime_error("build_from_edges: already built");
    if (edges_.empty()) throw std::runtime_error("build_from_edges: no edges appended");

    if (n_ <= 0)
    {
        throw std::invalid_argument("Graph::build_from_edges: n must be positive");
    }

    for (Edge edge : this->edges_)
    {
        printf("edge %d | u: %d v: %d w: %f\n", edge.id, edge.u, edge.v, edge.w);
    }

    std::vector<InputEdge> in;
    in.reserve(edges_.size());
    for (const auto& e : edges_)
        in.push_back({e.u, e.v, e.w});
    
    std::vector<Edge> canon;
    canonicaliseAndMerge(in, canon);
    if (canon.empty())
        throw std::runtime_error("build_from_edges: canonical edge list is empty.");
    
    
    std::vector<double> degree(n_, 0.0);
    for (const auto& e : canon)
    {
        degree[e.u] += e.w;
        degree[e.v] += e.w;
    }

    /* build CSR representation */
    std::vector<int>    row_ptr, col_idx;
    std::vector<double> adj_w;


    this->built_ = true;
}


void Graph::canonicaliseAndMerge (std::vector<InputEdge>& in, std::vector<Edge>& out)
{   
    std::cout << "in cononicaliseAndMerge\n";
    for (auto& e : in)
    {
        if (e.u < 0 || e.u >= n_ || e.v < 0 || e.v >= n_)
            throw std::invalid_argument("Edge endpoint out of range.");
        if (e.u == e.v)
            throw std::invalid_argument("Self-loop not allowed.");
        if (e.w <= 0.0)
            throw std::invalid_argument("Edge weight must be positive.");
        if (e.v < e.u)
            std::swap(e.u, e.v);
    }

    std::sort(in.begin(), in.end(), 
            [](const InputEdge& a, const InputEdge& b) /* sorting criterion */
            {
                return (a.u < b.u) || (a.u == b.u && a.v < b.v);
            }
    );

    out.clear();
    out.reserve(in.size());

    for (size_t i = 0; i < in.size();)
    {
        int    u = in[i].u;
        int    v = in[i].v;
        double wsum = 0.0;

        size_t j = i;
        while (j < in.size() && in[j].u == u && in[j].v == v)
        {
            wsum += in[j].w;
            ++j;
        }
        out.push_back({u, v, wsum, static_cast<int>(out.size())});
        i = j;
    }
}