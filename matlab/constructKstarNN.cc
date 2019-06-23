#include "armaMex.hpp"
#include "NetDECODE.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i;
	
	
	mat logPvals;
	int m = (int ) mxGetM(prhs[0]);
	int n = (int ) mxGetN(prhs[0]);
	if(mxIsSparse(prhs[0])) {
		logPvals = mat(armaGetSparseMatrix(prhs[0]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[0])); 
		
		logPvals = mat(ptr, m, n, true, true);		
	}

	double LC = 1.0, pval_threshold = 0.01;
	if(nrhs > 1) {
		LC = mxGetScalar(prhs[1]);
	}	
	if(nrhs > 2) {
		pval_threshold = mxGetScalar(prhs[2]);
	}			
	
	mat G = NetDECODEns::constructKstarNN(logPvals, LC, pval_threshold);


	plhs[0] = mxCreateDoubleMatrix(G.n_rows, G.n_cols, mxREAL);
	memcpy(mxGetPr(plhs[0]), G.memptr(), G.n_elem * sizeof(double)); 	
}
