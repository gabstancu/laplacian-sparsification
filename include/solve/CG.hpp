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
            
        }
};

#endif // CG_HPP