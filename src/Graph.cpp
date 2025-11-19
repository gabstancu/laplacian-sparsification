#include "graph/Graph.hpp"
#include <iostream>
#include "utils/display.hpp"

Graph::Graph(int n,
             std::vector<Edge>&& edges,
             std::vector<int>&& row_ptr,
             std::vector<int>&& col_idx,
             std::vector<double>&& adj_w,
             std::vector<double>&& degree,
             bool connected)
  : n_(n),
    m_(static_cast<int>(edges.size())),
    edges_(std::move(edges)),
    row_ptr_(std::move(row_ptr)),
    col_idx_(std::move(col_idx)),
    adj_w_(std::move(adj_w)),
    degree_(std::move(degree)),
    connected_(connected) {}



void Graph::build_from_edges ()
{   
    if (built_)         
        throw std::runtime_error("build_from_edges: already built");

    if (edges_.empty()) 
        throw std::runtime_error("build_from_edges: no edges appended");

    if (n_ <= 0)
    {
        throw std::invalid_argument("Graph::build_from_edges: n must be positive. Set n before building.");
    }

    std::vector<InputEdge> in;
    in.reserve(edges_.size());
    for (const auto& e : edges_)
    {   
        in.push_back({e.u, e.v, e.w});
    }
    
    std::vector<Edge> canon;
    canonicaliseAndMerge(in, canon);
    if (canon.empty())
    {
        throw std::runtime_error("build_from_edges: canonical edge list is empty.");
    }
    
    
    std::vector<double> degree(n_, 0.0);
    for (const auto& e : canon)
    {
        degree[e.u] += e.w;
        degree[e.v] += e.w;
    }

    /* build CSR representation */
    std::vector<int>    row_ptr, col_idx;
    std::vector<double> adj_w;
    buildCSR(canon, row_ptr, col_idx, adj_w);

    this->built_     = true;
    this->edges_     = std::move(canon);
    this->row_ptr_   = std::move(row_ptr);
    this->col_idx_   = std::move(col_idx);
    this->adj_w_     = std::move(adj_w);
    this->degree_    = std::move(degree);

    bool conn = check_connectivity();
    this->connected_ = conn;
}


void Graph::build_from_vertices ()
{
    if (built_)         
        throw std::runtime_error("build_from_vertices: already built");
    if (vertices_.empty()) throw std::runtime_error("build_from_vertices: no vertices appended");

    if (n_ <= 0)
    {
        throw std::invalid_argument("Graph::build_from_vertices: n must be positive. Set n before building.");
    }

    std::vector<InputEdge> in;
    for (const auto& vertex : vertices_)
    {
        for (const auto& [v, w] : vertex.neighbors)
        {
            in.push_back({vertex.u, v, w});
        }
    }

    std::vector<Edge> canon;
    canonicaliseAndMerge(in, canon);
    if (canon.empty())
    {
        throw std::runtime_error("build_from_vertices: canonical edge list is empty.");
    }

    std::vector<double> degree(n_, 0.0);
    for (const auto& e : canon)
    {
        degree[e.u] += e.w;
        degree[e.v] += e.w;
    }

    /* build CSR representation */
    std::vector<int>    row_ptr, col_idx;
    std::vector<double> adj_w;
    buildCSR(canon, row_ptr, col_idx, adj_w);

    this->built_     = true;
    this->edges_     = std::move(canon);
    this->row_ptr_   = std::move(row_ptr);
    this->col_idx_   = std::move(col_idx);
    this->adj_w_     = std::move(adj_w);
    this->degree_    = std::move(degree);

    bool conn = check_connectivity();
    this->connected_ = conn;
}



void Graph::buildCSR (std::vector<Edge>& edges, 
                      std::vector<int>& row_ptr, 
                      std::vector<int>& col_idx, 
                      std::vector<double>& adj_w)
{
    row_ptr.assign(n_+1, 0); /* cumulative offsets */

    for (const auto& e: edges) /* first using it to count neighbors */
    {
        // std::cout << "edge " << e.id << " e.u: " << e.u << " e.v: " << e.v << " weight: " << e.w << '\n';
        row_ptr[e.u + 1] += 1;
        row_ptr[e.v + 1] += 1;
    }

    // for (int a : row_ptr)
    // {
    //     std::cout << a << " ";
    // }
    // std::cout << '\n';

    for (int i = 0; i < n_; i++)
    {
        row_ptr[i + 1] += row_ptr[i];
    }

    // for (int a : row_ptr)
    // {
    //     std::cout << a << " ";
    // }
    // std::cout << '\n';

    col_idx.resize(row_ptr.back());
    adj_w.resize(row_ptr.back());

    std::vector<int> cursor = row_ptr;
    vertices_.resize(n_);
    std::vector<int> filled(n_, 0);

    for (const auto& e : edges)
    {
        int pos_u      = cursor[e.u];  // reserve next slot in row u
        col_idx[pos_u] = e.v;
        adj_w[pos_u]   = e.w;
        cursor[e.u]++;

        int pos_v      = cursor[e.v];  // symmetric entry in row v
        col_idx[pos_v] = e.u;
        adj_w[pos_v]   = e.w;
        cursor[e.v]++;

        if (filled[e.u])
        {
            vertices_[e.u].neighbors[e.v] = e.w;
        }
        else
        {
            vertices_[e.u] = {e.u, 0.0, {}}; filled[e.u] = 1;
        }

        if (filled[e.v])
        {
            vertices_[e.v].neighbors[e.u] = e.w;
        }
        else
        {
            vertices_[e.v] = {e.v, 0.0, {}}; filled[e.v] = 1;
        }
    }

    assert(static_cast<size_t>(row_ptr.back()) == 2 * edges.size());

    for (int u = 0; u < n_; ++u) 
    {
        assert(row_ptr[u] <= row_ptr[u + 1]);
    }

}



bool Graph::check_connectivity ()
{   
    if (n_ == 0)
    {
        return true;
    }

    std::vector<int> vis(n_, 0);
    std::queue<int>  q; 

    vis[0] = 1;
    q.push(0);
    int visited = 1;

    while (!q.empty())
    {
        /* for every node in q mark all it's neighbors as visited */
        int u = q.front(); 
        q.pop();
        for (int k = row_ptr_[u]; k < row_ptr_[u+1]; k++)
        {
            int v = col_idx_[k];
            if (!vis[v])
            {
                vis[v] = 1;
                visited++;
                q.push(v);
            }
        }        
    }

    return visited == n_;
}



void Graph::canonicaliseAndMerge (std::vector<InputEdge>& in, std::vector<Edge>& out)
{   
    // std::cout << "in cononicaliseAndMerge\n";
    for (auto& e : in)
    {
        if (e.u < 0 || e.u >= n_ || e.v < 0 || e.v >= n_)
            throw std::invalid_argument("Edge endpoint out of range.");
        if (e.u == e.v)
            throw std::invalid_argument("Self-loop not allowed.");
        if (e.w <= 0.0)
            throw std::invalid_argument("Edge weight must be positive.");
        if (e.v < e.u)
            std::swap(e.u, e.v);
    }

    std::sort(in.begin(), in.end(), 
            [](const InputEdge& a, const InputEdge& b) /* sorting criterion */
            {
                return (a.u < b.u) || (a.u == b.u && a.v < b.v);
            }
    );

    out.clear();
    out.reserve(in.size());

    for (size_t i = 0; i < in.size();)
    {
        int    u = in[i].u;
        int    v = in[i].v;
        double wsum = 0.0;

        size_t j = i;
        while (j < in.size() && in[j].u == u && in[j].v == v)
        {
            wsum += in[j].w;
            ++j;
        }
        out.push_back({u, v, wsum, static_cast<int>(out.size())});
        i = j;
    }
}


void Graph::print ()
{
    /* edges */

    /* nodes */

    /* general info */
}



SparseMatrix Graph::build_adjacency ()
{
    SparseMatrix A(n_, n_);
    std::vector<Triplet> EDGES;
    // EDGES.reserve(2 * m_);

    for (Edge e : edges_)
    {
        EDGES.emplace_back(e.u, e.v, e.w);
        EDGES.emplace_back(e.v, e.u, e.w);
    }

    A.setFromTriplets(EDGES.begin(), EDGES.end());
    A.makeCompressed();

    return A;
}



SparseMatrix Graph::buildLaplacianUnpinned ()
{
    SparseMatrix L(n_, n_);
    std::vector<Triplet> Entries;

    for (int i = 0; i < n_; i++) /* insert diagonal elements */
    {
        Entries.emplace_back(i, i, degree_[i]);
    }

    for (const auto& edge : edges_) /* insert off-diagonal elements */
    {
        Entries.emplace_back(edge.u, edge.v, -edge.w);
        Entries.emplace_back(edge.v, edge.u, -edge.w);
    }

    L.setFromTriplets(Entries.begin(), Entries.end());
    L.makeCompressed();

    return L;
}



SparseMatrix Graph::buildLaplacianPinned (int pinned_node, PinMaps* maps)
{
    SparseMatrix Lp(n_ - 1, n_ - 1);
    std::vector<Triplet> Entries;
    
    maps->pin.assign(n_, -1);
    
    for (int u = 0; u < n_; u++)
    {
        if (u == pinned_node)
            continue;

        maps->pin[u] = (u > pinned_node) ? maps->pin[u] = u - 1 : maps->pin[u] = u;
        maps->unpin.push_back(u);
    }

    // printVector(maps->pin, true);
    // printVector(maps->unpin, true);   

    for (int i = 0;  i < n_; i++) /* diagonal elements: keep their degree */
    {
        if (i == pinned_node)
            continue;

        Entries.emplace_back(maps->pin[i], maps->pin[i], degree_[i]);
    }

    for (const auto& edge: edges_)
    {
        if (edge.u == pinned_node || edge.v == pinned_node)
            continue;

        Entries.emplace_back(maps->pin[edge.u], maps->pin[edge.v], -edge.w);
        Entries.emplace_back(maps->pin[edge.v], maps->pin[edge.u], -edge.w);
    }

    Lp.setFromTriplets(Entries.begin(), Entries.end());
    Lp.makeCompressed();

    return Lp;
}