#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

#include "Eigen/Eigen"
#include <random>
#include <iostream>

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Vector       = Eigen::VectorXd;
class LinearSystem
{
    public:
        LinearSystem (SparseMatrix A, Vector b);
        LinearSystem () = default;

        void solve                        (); /* solver directly using Eigen's Sparse LU*/
        void estimate_extreme_eigenvalues (int max_iters = 1000, unsigned seed = 0);

        template <class Solver>
        void solve_with_solver (Solver& solver)
        {
            solver.solve(*this);
        }

        template <class Solver, class Preconditioner>
        void solve_with_solver (Solver& solver, Preconditioner& preconditioner)
        {
            solver.solve(*this, preconditioner);
        }

        void print_info ();

        /* get */
        SparseMatrix& A     () { return this->A_; }
        Vector&       b     () { return this->b_; }
        Vector&       u     () { return this->u_; }
        double        l_min () { return this-> lambda_min; }
        double        l_max () { return this-> lambda_max; }
        double        kappa () { return this-> kappa_hat; }
        int           n     () { return this->n_; }
        int           m     () { return this->m_; }

    private:
        int n_ = 0; // number of variables
        int m_ = 0; // number of equations

        double lambda_min = 0.0;
        double lambda_max = 0.0;
        double kappa_hat  = 0.0;

        SparseMatrix A_; // system matrix
        Vector       b_; // rhs vector
        Vector       u_; // solution vector

};

#endif // LINEAR_SYSTEM_HPP