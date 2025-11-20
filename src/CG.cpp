#include "solve/CG.hpp"

ConjugateGradient::ConjugateGradient (double tolerance, int max_iters)
{
    this->tol       = tolerance;
    this->max_iters = max_iters;
}


void ConjugateGradient::print ()
{
    
}