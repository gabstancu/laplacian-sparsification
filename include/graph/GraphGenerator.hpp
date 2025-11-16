#ifndef GRAPH_GENERATOR_HPP
#define GRAPH_GENERATOR_HPP

#include "graph/Graph.hpp"
#include <cstdint>
#include <random>

class GraphGenerator
{   

    public:
        struct Weights
        {
            enum class Kind 
            {
                Unit, Uniform, Lognormal, PowerLaw
            };
            Kind kind = Kind::Unit;

            double a     = 1.0, b     = 1.0; // uniform
            double mu    = 0.0, sigma = 1.0; // lognormal
            double alpha = 2.5, xmin  = 1.0; // power-law 
        };
        
        GraphGenerator () = default;

        /* Deterministic */
        Graph cycle (int n);
        Graph path  (int n);
        Graph grid  (int grid_size);
        Graph torus (int grid_size);
        /* Random */
        Graph random_d_regular (int n, int    d, bool ensure_connected);
        Graph erdos_renyi_gnp  (int n, double p, bool ensure_connected);
        Graph erdos_renyi_gnm  (int n, int    m, bool ensure_connected);

        void set_seed    (uint64_t s) { this->seed_ = s; };
        void set_weights (Weights w)  { this->w_    = w; }

    private:
        mutable int64_t seed_;
        std::mt19937_64 rng_ {0xC0FFEE};

        Weights w_;

        double add_weight_ ();

};


#endif // GRAPH_GENERATOR_HPP