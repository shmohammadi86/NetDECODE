#include <NetDECODE.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat asssessCoactivity(sp_mat &A) {
	
	mat logPvals = NetDECODEns::asssessCoactivity(A);
		
	return logPvals;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat constructKstarNN(mat &logPvals, double L_C = 5.0, double pval_threshold = 0.01) {

	mat G = NetDECODEns::constructKstarNN(logPvals, L_C, pval_threshold);
	
	return(G);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat symmetrizeNetwork(mat &G) {

	mat G_sym = NetDECODEns::symmetrizeNetwork(G);
	
	return(G_sym);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat predictActivityScores(sp_mat &A, IntegerVector rows, IntegerVector columns, mat &G) {	
	uvec rows_uvec(rows.size());
	for (int i = 0; i < rows_uvec.n_elem; i++) {
		rows_uvec(i) = (uword)(rows(i)-1);
	}	
	
	uvec cols_uvec(columns.size());
	for (int i = 0; i < cols_uvec.n_elem; i++) {
		cols_uvec(i) = (uword)(columns(i)-1);
	}		
	
	mat activity_profile = NetDECODEns::predictActivityScores(A, rows_uvec, cols_uvec, G);
	
	return(activity_profile);
}



