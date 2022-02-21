#ifndef STARTER_H
#define STARTER_H

//#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG

#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>

#include "armadillo"
using namespace arma;
using namespace std;

namespace NetDECODEns {
	mat asssessCoactivity(sp_mat &A);
	mat constructKstarNN(mat &logPvals, double L_C, double pval_threshold);
	mat symmetrizeNetwork(mat &G);
	mat predictActivityScores(sp_mat &A, uvec rows, uvec columns, mat &G);
}
#endif
