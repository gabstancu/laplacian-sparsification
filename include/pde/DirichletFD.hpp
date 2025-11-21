#ifndef DIRICHLET_FD_HPP
#define DIRICHLET_FD_HPP

#include "graph/Graph.hpp"

using Vector = Eigen::VectorXd;
namespace pde 
{   
    struct FDIndex 
    {
        int N = 0;             // grid is N x N, spacing h = 1/(N-1)
        inline int                id(int i, int j) const { return i * N + j; }
        inline std::pair<int,int> rc(int k)        const { return { k / N, k % N }; }
        inline double             h()              const { return (N > 1) ? 1.0 / (N - 1) : 1.0; }
    };

    struct DirichletBC
    {
        std::function<double(double)> left;   // g(0, y)
        std::function<double(double)> right;  // g(1, y)
        std::function<double(double)> bottom; // g(x, 0)
        std::function<double(double)> top;    // g(x, 1)
    };

    struct DirichletMaps
    {
        std::vector<int>  pin;
        std::vector<int>  unpin;
        std::vector<int>  is_boundary;
    };


    DirichletMaps build_dirichlet_maps (const FDIndex& index);

    
    SparseMatrix build_dirichlet_laplacian (Graph&         graph, 
                                            const DirichletMaps& maps); // build L_{II}
    
    // Assemble b_I = h^2 f(x,y) + Σ_{boundary nbr j} w_ij * g_j
    Vector build_dirichlet_rhs (Graph&         graph, 
                                const FDIndex&       idx, 
                                const DirichletMaps& maps,
                                const std::function<double(double,double)>& f, 
                                const DirichletBC&   g);

    // Lift interior solution back to full N^2 vector u, filling boundary by g
    Vector lift_dirichlet_solution (const FDIndex&       idx, 
                                    const DirichletMaps& maps, 
                                    const Vector&        u_I, 
                                    const DirichletBC&   g); // lift u_{I} to u
}


#endif // DIRICHLET_FD_HPP