#include <iostream>
#include "Eigen/Eigen"
#include "ginac.h"
#include "graph/Graph.hpp"
#include "graph/GraphGenerator.hpp"
#include "utils/display.hpp"
#include "pde/DirichletFD.hpp"
#include "LinearSystem.hpp"
#include "solve/CG.hpp"
#include <functional>

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
    Graph graph = generator.grid(5);

    pde::FDIndex       idx{5};
    pde::DirichletMaps maps      = pde::build_dirichlet_maps(idx);
    SparseMatrix       laplacian = pde::build_dirichlet_laplacian(graph, maps);


    pde::DirichletBC BC;
    BC.left   = [](double y) { return 0.0; };
    BC.right  = [](double y) { return 0.0; };
    BC.bottom = [](double x) { return std::sin(M_PI * x); };
    BC.top    = [](double x) { return 4.0 * std::sin(3.0 * M_PI * x); };

    std::function<double(double, double)> f = [](double x, double y)
    {
        return 0.0;
    };

    Vector b = pde::build_dirichlet_rhs(graph, idx, maps, f, BC); 

    LinearSystem system(laplacian, b);
    ConjugateGradient CG;
    system.solve_with_solver(CG);

    // printVector(graph.row_ptr(), true);
    // printVector(graph.col_idx(), true);
    // printVector(graph.adj_w(),   true);
    // printVector(graph.degree(),  true);


    return 0;
}