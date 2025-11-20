#include "pde/DirichletFD.hpp"

namespace pde 
{   
    using Triplet = Eigen::Triplet<double>;


    static inline bool is_boundary (int i, int j, int N)
    {
        return (i == 0) || (j == 0) || (i == N - 1) || (j == N -1);
    }


    static inline double boundary_node_value (int k, FDIndex& idx, DirichletBC& g)
    {
        int    N      = idx.N;
        auto   [i, j] = idx.rc(k);
        double h      = idx.h(); 
        double x      = j * h;
        double y      = 1.0 - i * h;

        if (i == 0)
            return g.top ? g.top(x) : 0.0;
        if (i == N - 1)
            return g.bottom ? g.bottom(x) : 0.0;
        if (j == 0)
            return g.left ? g.left(y) : 0.0;

        return g.right ? g.right(y) : 0.0;
    }


    DirichletMaps build_dirichlet_maps (FDIndex& index)
    {
        DirichletMaps maps;
        int N = index.N;
        int N2 = N * N;

        maps.is_boundary.resize(N2, 0);
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                int k = index.id(i, j);
                if (i == 0 || i == N - 1 || j == 0 || j == N - 1)
                {
                    maps.is_boundary[k] = 1;
                    maps.pin.push_back(k);
                }
                else
                {
                    maps.is_boundary[k] = 0;
                    maps.unpin.push_back(k);
                }
            }
        }

        return maps;
    }


    SparseMatrix build_dirichlet_laplacian (Graph&         graph, 
                                            DirichletMaps& maps)
    {
        const int n = static_cast<int>(maps.pin.size());

        const auto& row_ptr = graph.row_ptr();
        const auto& col_idx = graph.col_idx();
        const auto& adj_w   = graph.adj_w();
        const auto& degree  = graph.degree();  

        // Count interior nodes and roughly estimate nnz
        int nI = 0;
        for (int k = 0; k < n; ++k) 
            if (maps.pin[k] >= 0) ++nI;

        SparseMatrix         A(nI, nI);
        std::vector<Triplet> T;

        // 1 diagonal per interior + up to 4 off-diagonals on a grid
        T.reserve(nI * 5);

        for (int u = 0; u < n; ++u) 
        {
            const int ui = maps.pin[u];
            if (ui < 0) 
                continue; // boundary row skipped

            // Diagonal: full weighted degree (interior + boundary neighbors)
            T.emplace_back(ui, ui, degree[u]);

            // Off-diagonals: only interior neighbors contribute -w
            for (int k = row_ptr[u]; k < row_ptr[u + 1]; ++k) 
            {
                const int v  = col_idx[k];
                const int vi = maps.pin[v];
                if (vi >= 0) 
                {
                    T.emplace_back(ui, vi, -adj_w[k]);
                }
            }
        }

        A.setFromTriplets(T.begin(), T.end());
        A.makeCompressed();
        
        return A;
    }


    Vector build_dirichlet_rhs (Graph&         graph, 
                                FDIndex&       idx, 
                                DirichletMaps& maps,
                                std::function<double(double,double)>& f, 
                                DirichletBC&   g)
    {
        int    N = idx.N;
        int    n = static_cast<int>(maps.pin.size());
        double h = idx.h();

        const auto& row_prt = graph.row_ptr();
        const auto& col_idx = graph.col_idx();
        const auto& adj_w   = graph.adj_w();

        Vector b;
        b.setZero(static_cast<int>(maps.unpin.size()));

        for (int u = 0; u < n; u++)
        {
            int ui = maps.pin[u];

            if (ui < 0) /* boundary node */
            {
                continue;
            }

            const auto [i, j] = idx.rc(u);
            double x          = j * h;
            double y          = 1.0 - i * h;

            double rhs = 0.0;

            if (f)
                rhs -= (h*h) * f(x, y);

            for (int k = row_prt[u]; k < row_prt[u+1]; k++)
            {
                int v = col_idx[k];

                if (maps.pin[v] == -1)
                {
                    rhs += adj_w[v] * boundary_node_value(v, idx, g);
                }
            }
            b[ui] = rhs;
        }
        return b;
    }


    Vector lift_dirichlet_solution (FDIndex&       idx, 
                                    DirichletMaps& maps, 
                                    Vector&        u_I, 
                                    DirichletBC&   g)
    {   
        int N = idx.N;
        int n = N * N;

        Vector u(n);
        u.setZero();

        for (int k = 0; k < n; k++) /* fill boundaries */
        {
            if (maps.is_boundary[k])
            {
                u[k] = boundary_node_value(k, idx, g);
            } 
        }

        for (int ui = 0; ui < u_I.size(); ui++) /* fill interior */
        {
            int k = maps.unpin[ui];
            u[k]  = u_I[ui];
        }
        
        return u;
    }

} // namespace pde