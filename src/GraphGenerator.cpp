#include "graph/GraphGenerator.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>

double GraphGenerator::add_weight_ ()
{   
    rng_.seed(seed_++);
    
    switch (w_.kind)
    {
        case Weights::Kind::Unit :
            return 1.0;
        
        case Weights::Kind::Uniform: 
        {
            if (w_.a > w_.b)
                std::swap(w_.a, w_.b);

            std::uniform_real_distribution<double> U(w_.a, w_.b);
            double w = U(rng_);
            
            if (w <= 0.0)
                return std::nextafter(0.0, 1.0);
            else
                return w;
        }
        case Weights::Kind::Lognormal:
        {
            std::lognormal_distribution<double> LN(w_.mu, w_.sigma);
            double w = LN(rng_);

            if (w <= 0.0)
                return std::nextafter(0.0, 1.0);
            else
                return w;
        }
        case Weights::Kind::PowerLaw: 
        {
            std::uniform_real_distribution<double> U(0.0, 1.0);
            double u = std::min(1.0 - 1e-12, std::max(1e-12, U(rng_)));
            double w = w_.xmin / std::pow(1.0 - u, 1.0 / w_.alpha);

            return w;
        }
    }
    return 1.0;
}


Graph GraphGenerator::path (int n)
{
    Graph graph;
    graph.set_n(n);

    for (int i = 0; i < n - 1; i++)
    {
        // std::cout << "u: " << i << " v: " << i + 1 << " weight: " << add_weight_() << '\n';
        graph.add_edge({i, i+1, add_weight_()});
    }

    graph.build_from_edges();

    return graph;
}


Graph GraphGenerator::cycle (int n)
{
    Graph graph;
    graph.set_n(n);

    for (int i = 0; i < n; i++)
    {   
        if (i == n - 1) /* end connects to start */
        {   
            std::cout << "u: " << i << " v: " << 0 << " weight: " << add_weight_() << '\n';
            graph.add_edge({i, 0, add_weight_()});
        }
        else
        {   
            std::cout << "u: " << i << " v: " << i + 1 << " weight: " << add_weight_() << '\n';
            graph.add_edge({i, i+1, add_weight_()});
        }
    }

    // for (int i = 0; i < n; ++i) {
    //     int j = (i + 1) % n;
    //     // add each undirected edge once; Graph canonicalizes/merges
    //     if (i < j) 
    //     {   
    //         std::cout << "u: " << i << " v: " << j << " weight: " << add_weight_() << '\n';
    //         graph.add_edge({i, j, add_weight_()});
    //     }
        
    //     else
    //     {   
    //         std::cout << "u: " << j << " v: " << i << " weight: " << add_weight_() << '\n';
    //         graph.add_edge({j, i, add_weight_()});
    //     }       
    // }

    graph.build_from_edges();

    return graph;
}


Graph GraphGenerator::grid (int grid_size)
{
    if (grid_size < 2)
        throw std::invalid_argument("grid: grid_size must be >= 2");

    Graph graph;
    graph.set_n(int(grid_size*grid_size)); 

    for (int i = 0; i < grid_size; i++) 
    {
        for (int j = 0; j < grid_size; j++) 
        {
            int k = i * grid_size + j;
            if (j < grid_size - 1) 
            {
                graph.add_edge({k, k + 1, 1.0}); /* right */
            }
            if (i < grid_size - 1) 
            {
                graph.add_edge({k, k + grid_size, 1.0}); /* down */
            }
        }
    }

    graph.build_from_edges();
    return graph;
}


Graph GraphGenerator::torus (int grid_size)
{
    Graph graph;
    graph.set_n(int(grid_size*grid_size));

    // auto map = [grid_size](int i, int j)
    // {
    //     return i * grid_size + j;
    // };

    // auto fold = [grid_size](int a)
    // {

    // };

    for (int i = 0; i < grid_size; i++)
    {
        for (int j = 0; j < grid_size; j++)
        {
            int k      = i * grid_size + j;
            int right  = i * grid_size + (j + 1) % grid_size;
            int bottom = ((i+1) % grid_size) * grid_size + j;

            graph.add_edge({k, right, 1.0});
            graph.add_edge({k, bottom, 1.0});
        }
    }
    graph.build_from_edges();

    return graph;
}
