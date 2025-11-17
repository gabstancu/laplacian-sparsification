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
        LinearSystem (SparseMatrix A, Vector b) {};

        void solve                        (); /* solver directly using Eigen's Sparse LU*/
        void estimate_extreme_eigenvalues (int max_iters = 1000, unsigned seed);

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

    private:
        int n_ = 0; // number of variables
        int m_ = 0; // number of equations

        double lambda_min = 0.0;
        double lambda_max = 0.0;
        double kappa_hat  = 0.0;

        SparseMatrix    A_; // system matrix
        Eigen::VectorXd b_; // rhs vector
        Eigen::VectorXd u_; // solution vector

};

#endif // LINEAR_SYSTEM_HPP