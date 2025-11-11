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

        struct InputEdge /* appending */
        {
            int    u;
            int    v;
            double w;
        };

        struct PinMaps
        {
            std::vector<int> pin;
            std::vector<int> unpin;
        };
        
        static Graph FromEdges(int n, const std::vector<InputEdge>& in_edges);
        Graph () = default;

        void print ();

        int  n         () { return this->n_; }
        int  m         () { return this->m_; }
        bool connected () { return this->connected_; }

        void add_edge           (Edge edge) { this->edges_.push_back(edge); this->m_++;}
        void build_from_edges   ();
        std::vector<Edge> edges ()          { return this->edges_; }

        SparseMatrix build_adjacency        ();
        SparseMatrix build_degree           ();
        SparseMatrix buildLaplacianUnpinned ();
        SparseMatrix buildLaplacianPinned   ();


        std::vector<int>&    row_prt () { return this->row_ptr_; }
        std::vector<int>&    col_idx () { return this->col_idx_; }
        std::vector<double>& adj_w   () { return this->adj_w_; }

        std::pair<int*, int*>       neighborIndexRange  (int u) const;
        std::pair<double*, double*> neighbotWeightRange (int u) const;

        void set_n (int n) { this->n_ = n; }
        void set_m (int m) { this->m_ = m; }

    private:
        int n_ = 0; // nodes
        int m_ = 0; // edges

        std::vector<int>    row_ptr_;
        std::vector<int>    col_idx_;
        std::vector<double> adj_w_;
        std::vector<double> degree_; // cached

        std::vector<Edge> edges_;

        SparseMatrix A; // adjacency matrix
        SparseMatrix D; // degree    matrix

        bool connected_ = false;
        int  built_     = false;

        void buildCSR             ();
        bool check_connectivity   ();
        void canonicaliseAndMerge (std::vector<InputEdge>& in);

        Graph (int n, std::vector<Edge>&& edges, 
                      std::vector<int>&&  row_prt,
                      std::vector<int>&&  col_idx,
                      std::vector<double>&& adj_w,
                      std::vector<double>&& degree,
                      bool connected);
};

#endif // GRAPH_HPP