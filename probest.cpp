//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <omp.h>


// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
// R quantile
double quantilecpp(vec x,double q){
	vec tempq = sort(x);
	int sizex = tempq.n_elem;
	int indexq = round(sizex *q	);
	double out = tempq(indexq-1);
	return(out);
}

// [[Rcpp::export]]
List probest(const mat X, const mat beta, const int cores){
	
//Estimate mean and credible bounds of probability (logit)
	//Function arguments
		//X - covariate matrix
		//beta - parameter matrix, row = MCMC draw
		//cores - number of cores to use for parallell processing
		
		//container for probabilities
		mat prob(beta.n_rows, X.n_rows, fill::zeros);
		
		//calculate probabilities
		#pragma omp parallel for
		for(int i=0; i<beta.n_rows; ++i){
			for(int j=0; j<X.n_rows; ++j){
				double temp = as_scalar(exp(X(j,span::all) * beta(i,span::all).t()));
				prob(i,j) = temp /(1 + temp);
			}	
		}
		
		//summary stats
		
		mat means = mean(prob, 0); //posterior mean @ each covariate value
		vec lcl(prob.n_cols, fill::zeros);
		vec ucl(prob.n_cols, fill::zeros);
		
		for(int i=0; i<prob.n_cols; ++i){
			vec temp =  conv_to< colvec >::from(prob(span::all, i));
			lcl(i) = quantilecpp(temp, 0.025);
			ucl(i) = quantilecpp(temp, 0.975);
		}
		
		
		
		return List::create(Rcpp::Named("mu") = means,
							Rcpp::Named("quantile025") = lcl,
							Rcpp::Named("quantile975") = ucl
							);
		
		
}
