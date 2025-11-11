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

    // if (n_ <= 0)
    // {
    //     throw std::invalid_argument("Graph::build_from_edges(): n must be positive");
    // }

    for (Edge edge : this->edges_)
    {
        printf("edge %d | u: %d v: %d w: %f\n", edge.id, edge.u, edge.v, edge.w);
    }

    std::vector<InputEdge> in;
    in.reserve(edges_.size());
    for (const auto& e : edges_)
        in.push_back({e.u, e.v, e.w});
    
    canonicaliseAndMerge(in);
    
    


    this->built_ = true;
}


void Graph::canonicaliseAndMerge (std::vector<InputEdge>& in)
{   
    std::cout << "in conicaliseAndMerge\n";
    for (auto& e : in)
    {
        
    }
}