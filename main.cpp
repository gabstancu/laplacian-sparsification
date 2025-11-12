#include <iostream>
#include "Eigen/Eigen"
#include "graph/Graph.hpp"

int main ()
{   
    Graph::Edge edge_1{0, 1, 1.0, 0};
    Graph::Edge edge_2{0, 2, 1.0, 1};
    Graph::Edge edge_3{1, 2, 1.0, 2};
    Graph::Edge edge_4{2, 3, 1.0, 3};
    Graph::Edge edge_5{0, 3, 1.0, 4};

    Graph graph;
    graph.set_n(4);
    graph.add_edge(edge_1);
    graph.add_edge(edge_2);
    graph.add_edge(edge_3);
    graph.add_edge(edge_4);
    graph.add_edge(edge_5);
    graph.build_from_edges();

    return 0;
}