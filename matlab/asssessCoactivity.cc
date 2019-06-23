#include "armaMex.hpp"
#include "NetDECODE.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i;
	
	
	sp_mat A;
	int m = (int ) mxGetM(prhs[0]);
	int n = (int ) mxGetN(prhs[0]);
	if(mxIsSparse(prhs[0])) {
		A = armaGetSparseMatrix(prhs[0]);  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[0])); 
		
		A = sp_mat(mat(ptr, m, n, true, true));		
	}

	mat logPvals = NetDECODEns::asssessCoactivity(A);


	plhs[0] = mxCreateDoubleMatrix(logPvals.n_rows, logPvals.n_cols, mxREAL);
	memcpy(mxGetPr(plhs[0]), logPvals.memptr(), logPvals.n_elem * sizeof(double)); 	
}
