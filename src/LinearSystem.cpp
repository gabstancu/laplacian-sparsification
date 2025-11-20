#include "LinearSystem.hpp"

LinearSystem::LinearSystem (SparseMatrix A, Vector b) : A_(std::move(A)), b_(std::move(b))
{
    if (A_.rows() != b_.size())
    {
        throw std::invalid_argument("Matrix A and vector b dimensions do not match.");
    }
    
    u_.setZero(b_.size());
    n_ = static_cast<int>(A_.cols());
    m_ = static_cast<int>(A_.rows());

    estimate_extreme_eigenvalues();
}


void LinearSystem::solve ()
{
    Eigen::SparseLU<SparseMatrix> solver;
    solver.analyzePattern(A_);
    solver.factorize(A_);
    u_ = solver.solve(b_);
}


void LinearSystem::estimate_extreme_eigenvalues (int max_iters, unsigned seed)
{
    int n = A_.rows();
    if (n == 0)
    {
        throw std::runtime_error("Empty matrix in Lanczos.");
    }

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dist(-1.0, 1.0);

    Vector q = Vector::NullaryExpr(n, [&]() { return dist(gen); } );
    q.normalize();

    Vector q_prev = Vector::Zero(n);
    Vector Aq     = Vector::Zero(n);

    std::vector<double> alpha, beta;
    alpha.reserve(max_iters);
    beta.reserve(max_iters);

    for (int j = 0; j < max_iters; ++j)
    {
        Aq.noalias() = A_ * q;

        double a_j = q.dot(Aq);
        alpha.push_back(a_j);

        if (j > 0) 
            Aq.noalias() -= beta.back() * q_prev;
        Aq.noalias() -= a_j * q;

        double b_j = Aq.norm();
        if (b_j < 1e-14) break;

        beta.push_back(b_j);
        q_prev = q;
        q      = Aq / b_j;
    }

    // build the tridiagonal T
    const int m = static_cast<int>(alpha.size());
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m, m);

    for (int i = 0; i < m; ++i) 
    {
        T(i,i) = alpha[i];
        // only set off-diagonals if both indices are inside the matrix
        if (i+1 < m && i < static_cast<int>(beta.size())) 
        {
            T(i, i+1) = beta[i];
            T(i+1, i) = beta[i];
        }
    }

    if (m == 0) 
        throw std::runtime_error("Lanczos produced no alpha; check input matrix or max_steps.");

    if (m == 1) 
    {
        /* both extremes are the same */
        this->lambda_min = this->lambda_max = alpha[0];
        return;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(T);
    lambda_min = eigs.eigenvalues()(0);
    lambda_max = eigs.eigenvalues()(m-1);
    kappa_hat  = this->lambda_max / this->lambda_min; 
}


void LinearSystem::calc_omega_()
{
    this->omega_ = 2.0 / (1 + std::sin((M_PI) / this->n_ ));
    this->omega_ = std::min(this->omega_, 1.9 + (0.05) / this->n_);
}


void LinearSystem::print_info ()
{
    std::cout << "Linear System Info:\n";
    std::cout << "Number of variables (n): " << n_ << '\n';
    std::cout << "Number of equations (m): " << m_ << '\n';
    std::cout << "Estimated min eigenvalue: " << lambda_min << '\n';
    std::cout << "Estimated max eigenvalue: " << lambda_max << '\n';
    std::cout << "Estimated condition number (kappa_hat): " << kappa_hat << '\n';
}