#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "Eigen/Eigen"
#include "tuple"

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Triplet      = Eigen::Triplet<double>;


class Graph
{   
    public:
        struct Edge /* canonical edge (stored): u < v */
        {
            int    u;
            int    v;
            double w;
            int    id;
        };

        struct Vertex
        {
            int    u;
            double value;
            std::unordered_map<int, double> neighbors; // v -> w
        };

        struct InputEdge /* appending */
        {
            int    u;
            int    v;
            double w;
        };

        struct PinMaps
        {
            std::vector<int> pin;   /* -1 if node v is pinned */
            std::vector<int> unpin; /* remove element v-1 */
        };
        
        Graph () = default;

        void print ();

        int  n           () { return this->n_; }
        int  m           () { return this->m_; }
        bool connected   () { return this->connected_; }
        std::string type () { return this->type_; }

        void              add_edge            (Edge edge)     { this->edges_.push_back(edge);      this->m_++;}
        void              add_vertex          (Vertex vertex) { this->vertices_.push_back(vertex); this->n_++;}   
        void              build_from_edges    ();
        void              build_from_vertices ();
        std::vector<Edge> edges               ()          { return this->edges_; }
        std::vector<Vertex> vertices          ()          { return this->vertices_; }

        SparseMatrix build_adjacency        ();
        SparseMatrix buildLaplacianUnpinned ();
        SparseMatrix buildLaplacianPinned   (int pinned_node, PinMaps* maps);

        std::vector<int>&    row_ptr () { return this->row_ptr_; }
        std::vector<int>&    col_idx () { return this->col_idx_; }
        std::vector<double>& adj_w   () { return this->adj_w_; }
        std::vector<double>& degree  () { return this->degree_; }

        void set_n    (int n)                   { this->n_ = n; }
        void set_m    (int m)                   { this->m_ = m; }
        void set_type (const std::string& type) { this->type_ = type; }

    private:
        int n_ = 0; // nodes
        int m_ = 0; // edges

        std::string type_ = "";

        std::vector<int>    row_ptr_;
        std::vector<int>    col_idx_;
        std::vector<double> adj_w_;
        std::vector<double> degree_; // cached

        std::vector<Edge>   edges_;
        std::vector<Vertex> vertices_;

        SparseMatrix A; // adjacency matrix
        SparseMatrix D; // degree    matrix

        bool connected_ = false;
        int  built_     = false;

        void buildCSR             (std::vector<Edge>& edges, 
                                   std::vector<int>& row_prt, 
                                   std::vector<int>& col_idx, 
                                   std::vector<double>& adj_w);
        bool check_connectivity   ();
        void canonicaliseAndMerge (std::vector<InputEdge>& in, std::vector<Edge>& out);

        Graph (int n, std::vector<Edge>&& edges, 
                      std::vector<int>&&  row_ptr,
                      std::vector<int>&&  col_idx,
                      std::vector<double>&& adj_w,
                      std::vector<double>&& degree,
                      bool connected);
};

#endif // GRAPH_HPP