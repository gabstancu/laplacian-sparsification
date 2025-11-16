#include "graph/GraphGenerator.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <numeric>
#include <random>

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
                graph.add_edge({k, k + 1, add_weight_()}); /* right */
            }
            if (i < grid_size - 1) 
            {
                graph.add_edge({k, k + grid_size, add_weight_()}); /* down */
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

            graph.add_edge({k, right, add_weight_()});
            graph.add_edge({k, bottom, add_weight_()});
        }
    }
    graph.build_from_edges();

    return graph;
}


Graph GraphGenerator::erdos_renyi_gnp (int n, double p, bool ensure_connected)
{   
    if (n < 1)  
        throw std::invalid_argument("Gnp: n must be >= 1");

    if (p < 0.0 || p > 1.0) 
        throw std::invalid_argument("Gnp: p must be in [0,1]");

    int max_tries = ensure_connected ? 12 : 1;   

    for (int t = 0; t < max_tries; t++)
    {   
        Graph graph;
        graph.set_n(n);

        rng_.seed(seed_++);
        std::bernoulli_distribution flip(p);

        for (int u = 0; u < n; u++) /* add each u < v with prob p */
        {
            for (int v = u + 1; v < n; v++)
            {
                // std::cout << "u: " << u << " v: " << v << '\n';
                int f = flip(rng_);
                if (f)
                {
                    // std::cout << "+ u: " << u << " v: " << v << '\n';
                    // std::cout << "flip(rng_): " << f << '\n';
                    graph.add_edge({u, v, add_weight_()});
                }
                // else
                // {   
                //     std::cout << "- u: " << u << " v: " << v << '\n';
                //     std::cout << "flip(rng_): " << f << '\n';
                // }
                // std::cout << "----------------------------\n";
            }
        }

        graph.build_from_edges();

        if (!ensure_connected || graph.connected())
        {
            return graph;
        }
    }
    throw std::runtime_error("Gnm: failed to obtain a connected graph within attempts; adjust n (or p) or disable ensure_connected.");
}


Graph GraphGenerator::erdos_renyi_gnm (int n, int m, bool ensure_connected)
{   
    if (n < 1)
        throw std::invalid_argument("Gnm: n must be >= 1");
    
    if (m < 0 || m > 1LL * n * (n - 1) / 2)
        throw std::invalid_argument("Gnm: m out of feasible range");

    auto key = [n](int u, int v) -> uint64_t /* hashset entry */
    {
        if (u > v) 
            std::swap(u, v);     
        return (uint64_t)u * (uint64_t)n + (uint64_t)v;
    };
    
    int max_tries = ensure_connected ? 12 : 1;

    for (int t = 0; t < max_tries; t++)
    {
        rng_.seed(seed_++);

        std::uniform_int_distribution<int> U(0, n-1);
        std::unordered_set<uint64_t> chosen;
        chosen.reserve((size_t)m * 2);

        while ((int)chosen.size() < m)
        {
            int u = U(rng_);
            int v = U(rng_);
            if (u == v)
                continue;
            chosen.insert(key(u, v));
        }

        Graph graph;
        graph.set_n(n);

        for (auto k : chosen)
        {
            int u = (int)(k / (uint64_t)n);
            int v = (int)(k % (uint64_t)n);
            graph.add_edge({u, v, add_weight_()});
        }

        graph.build_from_edges();

        if (!ensure_connected || graph.connected())
            return graph;
    }
    throw std::runtime_error("Gnm: failed to obtain a connected graph within attempts; adjust m or disable ensure_connected.");
}


Graph GraphGenerator::random_d_regular (int n, int d, bool ensure_connected)
{   
    if (n <= 0) 
        throw std::invalid_argument("d-regular: n must be > 0");
    if (d <= 0 || d >= n) 
        throw std::invalid_argument("d-regular: require 1 <= d < n");
    if ((n * 1LL * d) % 2 != 0) 
        throw std::invalid_argument("d-regular: n*d must be even");

    int max_tries = 50;

    for (int t = 0; t < max_tries; t++)
    {
        rng_.seed(seed_++);

        /* each vertex appears d times */
        std::vector<int> stubs;
        stubs.reserve(n*d);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < d; j++)
            {
                stubs.push_back(i);
            }
        }
        std::shuffle(stubs.begin(), stubs.end(), rng_);

        bool bad = false;
        Graph graph;
        graph.set_n(n);

        for (size_t k = 0; k < stubs.size(); k+=2)
        {
            int u = stubs[k], v = stubs[k+1];
            if (u == v)
            {
                bad = true;
                break;
            }
            graph.add_edge({u, v, add_weight_()});
        }

        if (bad)
            continue;

        graph.build_from_edges();

        auto& row_ptr = graph.row_ptr();
        bool ok = true;

        for (int i = 0; i < n; i++)
        {
            if (row_ptr[i+1] - row_ptr[i] != d)
            {
                ok = false;
                break;
            }
        }

        if (!ok)
            continue;

        if (ensure_connected && !graph.connected()) 
            continue;
        
        return graph;
    }  
    throw std::runtime_error("d-regular: failed to sample a simple (and connected) graph after many attempts; try smaller d or n.");
}