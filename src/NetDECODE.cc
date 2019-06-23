#include <NetDECODE.h>

namespace NetDECODEns {
	mat asssessCoactivity(sp_mat &A) {
		int Nr = A.n_rows;
		int Nc = A.n_cols;
		printf("NetDECODE::\n");
		sp_mat B = spones(A);
		
		vec p_j = vec(trans(mean(B, 0)));
		vec p_i = vec(mean(B, 1));
		double p = sum(p_i)/Nr;
				
		
		printf("\tComputing observed statistics\n"); fflush(stdout);
		mat obs = mat(B * trans(B));

		printf("\tComputing expected statistics\n"); fflush(stdout);		
		double rho = (norm(p_j, 2) / p);
		mat expect = (rho*rho)* (p_i * trans(p_i));
	
		uvec idx;
		mat Delta = (obs / expect) - 1;
		
		idx = find(expect == 0);
		Delta(idx).zeros();
		
		mat logPvals = ( (Delta % Delta) / (2 + Delta) ) % expect;

		idx = find(Delta < 0);
		logPvals(idx).zeros();
		
		logPvals.diag().zeros();		
		
		logPvals /= log(10); // Change log-base
		
		return logPvals;
	}

	
	// From:: Anava, O. & Levy, K. k\ast -Nearest Neighbors: From Global to Local. in Advances in Neural Information Processing Systems 29 (eds. Lee, D. D., Sugiyama, M., Luxburg, U. V, Guyon, I. & Garnett, R.) 4916–4924 (Curran Associates, Inc., 2016).
	mat constructKstarNN(mat &logPvals, double L_C = 1.0, double pval_threshold = 0.01) {
		printf("Running k*-nearest neighbors algorithm\n");
		
		logPvals = clamp(logPvals, 0, 300);		
		mat D = exp(-logPvals); // p-values		
				
		mat G_prime = zeros(size(logPvals));
		
		int gene_no = logPvals.n_rows;
		for(int i = 0; i < gene_no; i++) {
			vec d = D.col(i);
			double counts = sum(d < 1);
			
			if(counts == 0) {
				continue;
			}
			
			uvec perm = sort_index(d, "ascend");			
			vec beta = L_C * d(perm);		
			
			double k = 0, Sum_beta = 0, Sum_beta_square = 0;			
			double lambda = beta(0)+1;
			
			double last_lambda = lambda;
			for(k = 0; k < counts; k++) {
				if(lambda <= beta(k))
					break;
				else {
					last_lambda = lambda;
				}
				
				Sum_beta += beta(k);
				Sum_beta_square += std::pow(beta(k), 2);

				lambda = (1.0 / k) * ( Sum_beta + sqrt( k  + std::pow(Sum_beta, 2) - k * Sum_beta_square ) ) ; 				
			}
			int knn = (int)(k-1);
						
			vec w = last_lambda-beta(span(0, knn));//last_lambda - beta(span(0, knn));
			w /= sum(w);
									
			vec v = zeros(gene_no);
			v(perm(span(0, knn))) = w;
			
			G_prime.col(i) = v;
		}
		
		return(G_prime);
	}
	
	// Uses geometric mean of the edge weights
	mat symmetrizeNetwork(mat &G) {
		mat Gt = trans(G);		
		mat G_sym = sqrt(G % Gt);		
		
		return(G_sym);
	}
	
	
	mat predictActivityScores(sp_mat &A, uvec rows, uvec cols, mat &G) {	
		mat logPvals = zeros(rows.n_elem, cols.n_elem);
		
		sp_mat B = spones(A);
		
		vec p_c = vec(trans(mean(B, 0)));

		mat Gt = trans(G.cols(rows));

		// Probabilistic Construction of Deterministic Algorithms: Approximating packing integer programs (1988), Section 1.1
		for(int j = 0; j < cols.n_elem; j++) {
			vec v = vec(B.col(cols(j)));
			
			double p = p_c(cols(j));
			
			vec Obs = Gt * v;
			vec Exp = p * sum(Gt, 1);

			vec Delta = (Obs / Exp) - 1;
			
			for(int i = 0; i < rows.n_elem; i++) {				
				logPvals(i, j) = Exp(i) * (pow(1 + Delta(i), 1+Delta(i)) - Delta(i));
			}
		}

		return(logPvals);
	}
	
}
