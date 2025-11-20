#ifndef CG_HPP
#define CG_HPP

#include "Eigen/Eigen"

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using DenseMatrix  = Eigen::MatrixXd;
using Vector       = Eigen::VectorXd;

class ConjugateGradient
{
    private:
        double      tol       = 1e-6;
        int         max_iters = 1e6;
        std::string name      = "CG";
        Vector      u_;

    public:
        ConjugateGradient () = default;
        ConjugateGradient (double tolerance, int max_iterations);

        template<class System>
        void solve (System& system)
        {
            const auto& A = system.A();
            const auto& b = system.b();
                  auto& u = system.u();


            Vector r      = b - A * u;
            double b_norm = b.norm();
            double r_norm = r.norm();

            if (r_norm / b_norm <= tol)
            {
                u_ = u;
                return;
            }

            Vector d = r;
            Vector Ad(A.rows());

            for (int i = 0; i < max_iters; i++)
            {   
                // std::cout << "--------------------- iter. " << i+1 << " ---------------------\n";
                Ad.noalias() = A * d;

                double alpha      = ((r.transpose() * r) / (d.transpose() * Ad)).coeff(0); // step size
                double r_prev_dot = (r.transpose() * r).coeff(0); // to calculate beta

                u.noalias() += alpha * d;
                r.noalias() -= alpha * Ad;

                r_norm = r.norm();

                if (r_norm / b_norm <= tol)
                {
                    u_ = u;
                    return;
                }

                double beta = r.dot(r) / r_prev_dot;
                d           = r + beta * d; // update direction
            }
            u_ = u;
            return;
        }

        void print ();

        void set_tolerance      (double tolerance)   { this->tol       = tolerance; }
        void set_max_iterations (int max_iterations) { this->max_iters = max_iterations; }

        double tolerance      () { return this->tol;}
        int    max_iterations () { return this-> max_iters; }
};

#endif // CG_HPP