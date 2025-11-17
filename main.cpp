#include <iostream>
#include "Eigen/Eigen"
#include "graph/Graph.hpp"
#include "graph/GraphGenerator.hpp"
#include "utils/display.hpp"

int main ()
{
    // Graph::Edge edge_1{0, 1, 1.0};
    // Graph::Edge edge_2{0, 2, 1.0};
    // Graph::Edge edge_3{1, 2, 1.0};
    // Graph::Edge edge_4{2, 3, 1.0};
    // Graph::Edge edge_5{0, 3, 1.0};

    // Graph graph;
    // graph.set_n(4);
    // graph.add_edge(edge_1);
    // graph.add_edge(edge_2);
    // graph.add_edge(edge_3);
    // graph.add_edge(edge_4);
    // graph.add_edge(edge_5);
    // graph.build_from_edges();

    // Graph::PinMaps maps;
    // std::cout << graph.buildLaplacianPinned(0, &maps) << '\n';

    GraphGenerator generator;
    Graph graph = generator.grid(3);
    graph.build_adjacency();
    // SparseMatrix L = graph.buildLaplacianUnpinned();

    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << "i: " << i << ", j: " << j << '\n';
    //         std::cout << "k: " << i * 4 + j << '\n';

    //         if (j == 0) /* left */
    //         {
    //             std::cout << "left edge\n";
    //         }
    //         if (i == 0) /* top */
    //         {
    //             std::cout << "top edge\n";
    //         }
            
    //         if (j == 3) /* right */
    //         {
    //             std::cout << "right edge\n";
    //         }

    //         if (i == 3) /* bottom */
    //         {
    //             std::cout << "bottom edge\n";
    //         }
            
    //     }
    // }

    printVector(graph.row_ptr(), true);
    printVector(graph.col_idx(), true);
    printVector(graph.adj_w(), true);
    printVector(graph.degree(), true);


    return 0;
}