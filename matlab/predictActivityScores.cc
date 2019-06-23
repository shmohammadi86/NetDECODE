#include "armaMex.hpp"
#include "NetDECODE.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i;
	
	int vector_size;
	
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
	
	m = (int ) mxGetM(prhs[1]);
	n = (int ) mxGetN(prhs[1]);	
	if(min(m, n) != 1) {
		mexErrMsgTxt("Cols parameter should be a one dimensional vector.");
	}	
	vector_size = max(m, n);
	uvec rows(vector_size);
	if(mxIsSparse(prhs[1])) {
		mat temp = mat(armaGetSparseMatrix(prhs[1]));  
		for (i = 0; i < vector_size; i++) {
			rows(i) = (uword)(temp[i]-1);
		}		
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[1])); 
		for (i = 0; i < vector_size; i++) {
			rows(i) = (uword)(ptr[i]-1);
		}		
	}
	
	m = (int ) mxGetM(prhs[2]);
	n = (int ) mxGetN(prhs[2]);	
	if(min(m, n) != 1) {
		mexErrMsgTxt("Null_cols parameter should be a one dimensional vector.");
	}	
	vector_size = max(m, n);
	uvec cols(vector_size);
	if(mxIsSparse(prhs[1])) {
		mat temp = mat(armaGetSparseMatrix(prhs[1]));  
		for (i = 0; i < vector_size; i++) {
			cols(i) = (uword)(temp[i]-1);
		}		
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[1])); 
		for (i = 0; i < vector_size; i++) {
			cols(i) = (uword)(ptr[i]-1);
		}		
	}	
	
	
	mat Net;
	m = (int ) mxGetM(prhs[3]);
	n = (int ) mxGetN(prhs[3]);
	if(mxIsSparse(prhs[3])) {
		Net = mat(armaGetSparseMatrix(prhs[3]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[3])); 
		
		Net = mat(ptr, m, n, true, true);		
	}	

	mat activity_profile = NetDECODEns::predictActivityScores(A, rows, cols, Net);

	plhs[0] = mxCreateDoubleMatrix(activity_profile.n_rows, activity_profile.n_cols, mxREAL);
	memcpy(mxGetPr(plhs[0]), activity_profile.memptr(), activity_profile.n_elem * sizeof(double)); 		
}
