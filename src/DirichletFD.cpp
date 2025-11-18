#include "pde/DirichletFD.hpp"

namespace pde 
{   
    using Triplet = Eigen::Triplet<double>;



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
        for (int k = 0; k < n; ++k) if (maps.pin[k] >= 0) ++nI;
        SparseMatrix A(nI, nI);
        std::vector<Triplet> T;

        // 1 diagonal per interior + up to 4 off-diagonals on a grid
        T.reserve(nI * 5);

        for (int u = 0; u < n; ++u) {
            const int ui = maps.pin[u];
            if (ui < 0) 
                continue; // boundary row skipped

            // Diagonal: full weighted degree (interior + boundary neighbors)
            T.emplace_back(ui, ui, degree[u]);

            // Off-diagonals: only interior neighbors contribute -w
            for (int k = row_ptr[u]; k < row_ptr[u + 1]; ++k) {
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

} // namespace pde