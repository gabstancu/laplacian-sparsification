#ifndef DSB_HPP
#define DSB_HPP

namespace sparsify
{
    struct DSBStats
    {
        int      n         = 0;
        int      m_in      = 0; /* unique undirected edges seen */
        int      m_out     = 0; /* edges kept after sampling */
        double   p_min     = 0.0;
        double   p_mean    = 0.0;
        double   p_max     = 0.0;
        double   w_sum_in  = 0.0;
        double   w_sum_out = 0.0;
        uint64_t seed      = 0;
    };

/**
 * Degree-based importance sampling (DSB).
 * - Scores: s_e = w_e * (1/deg(u) + 1/deg(v))  (deg is *weighted* degree)
 * - Choose α by bisection so that sum_e min(1, α s_e) \approx m_target (capped at |E|)
 * - Sample edges independently with p_e = min(1, α s_e)
 * - Reweight kept edge to w_e / p_e  (unbiased: E[\tilde L] = L)
 *
 * @param  G         Input (undirected) graph (CSR built)
 * @param  m_target  Target expected number of kept edges
 * @param  seed      RNG seed
 * @param  out       Optional stats
 * @return sparsified graph (reweighted)
 */
 Graph dsb (const Graph& G, int m_target, uint64_t seed, DSBStats* out = nullptr);
    
} // sparsify

#endif // DSB_HPP