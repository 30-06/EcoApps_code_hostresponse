//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h>



// [[Rcpp::plugins(openmp)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


//Necessary functions
//Find interval of a vector of break points containing a number - Hadley advanced R 
IntegerVector findInterval2(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());

  NumericVector::iterator x_it = x.begin(), x_end = x.end(),
    breaks_it = breaks.begin(), breaks_end = breaks.end();
  IntegerVector::iterator out_it = out.begin(), out_end = out.end();
  NumericVector::iterator ubound; 
  
  for(; x_it != x_end; x_it++, out_it++) {
    ubound = std::upper_bound(breaks_it, breaks_end, *x_it);
    *out_it = std::distance(ubound, breaks_it);
  }

  return out;
}

//Based on code from: https://github.com/RcppCore/RcppArmadillo/blob/master/inst/include/RcppArmadilloExtensions/sample.h

//[[Rcpp::export]]
 int ProbSampleReplace(vec prob){ 
			double t1 = sum(prob);
			if(t1 !=0){
					prob=prob/t1;
			}
			int out;
            double rU;
            int ii, jj;
            int nOrig = prob.n_elem;
            int nOrig_1 = nOrig - 1;
            uvec perm = sort_index(prob, "descend"); //descending sort of index
            prob = sort(prob, "descend");  // descending sort of prob
            
            // cumulative probabilities 
            prob = cumsum(prob);
            // compute the sample 
             rU = unif_rand();
                for (jj = 0; jj < nOrig_1; jj++) {
                    if (rU <= prob[jj])
                        break;
                }
            out = perm[jj];
			return(out);
 }


// [[Rcpp::export]]
// R quantile
double quantilecpp(vec x,double q){
	vec tempq = sort(x);
	int sizex = tempq.n_elem;
	int indexq = round(sizex *q	);
	double out = tempq(indexq-1);
	return(out);
}

//[[Rcpp::export]]
//determine failure day

double fail(const vec surv){
	int out;
	vec temp1(surv.n_elem+1, fill::zeros);
	temp1(0) = 1-surv(0);
	//Rcpp::Rcout<<"fail"<< std::endl;
	
	for(uword i=1; i<surv.n_elem; ++i){  //assumes start with day 0
		temp1(i) = as_scalar((surv(i-1)-surv(i))); //daily prob of failing
	}

	temp1(surv.n_elem) = 1-sum(temp1);
	out = ProbSampleReplace(temp1); //fail day
		//Rcpp::Rcout<<temp1<< std::endl;
	return(out);
}	



//Main Function

//[[Rcpp::export]] 

List deriveval(const vec gamma1, const vec gamma2, const vec gamma3, const vec gamma4,
		const vec gamma6, const mat beta01, const mat beta06, const mat beta0dd, 
		const mat beta0rec, const mat betah, const mat betah3, const mat betah4,
		const mat Ka, const mat Kadd, const int N, const vec days,const mat XX,const int cores){
			
		//gamma1 - baseline hazard for dying - not disease
		//gamma2 - baseline hazard for first infection
		//gamma3 - baseline hazard for disease death
		//gamma4 - baseline hazard for recovery
		//gamma6 - baseline hazard for reinfection
		
		//beta01 - probability of infection - not disease (logit scale)
		//beta06 - probability of reinfection (logit scale)
		//beta0dd - probability of disease death (logit scale)
		//beta0rec - probability of recovery (logit scale)
		//betah - hazard suscept to infected parms
		//betah3 - hazard infected to disease death parms
		//betah4 - hazard infected to recovered parms
		//Ka - kernel time effect for infection hazard
		//Kadd - kernel time effect for disease death hazard
		
		//N - population size to simulate
		//days - day to estimate prop of N in each state
		//XX - design matrix
		//cores - number of cores
		
		omp_set_num_threads(cores);
		
		//Create containers for output
		
		mat Ssucept(Ka.n_rows, Ka.n_cols, fill::zeros); //survival curves for susceptible state
		mat Sinfect(Kadd.n_rows, Kadd.n_cols, fill::zeros); //survival curves for infected state
		mat Srec(Kadd.n_rows, Kadd.n_cols, fill::zeros); //survival curves for recovered state
		
		mat hazinf(Ka.n_rows, Ka.n_cols, fill::zeros); //hazard curve for susceptible to infected state
		mat hazdd(Kadd.n_rows, Kadd.n_cols, fill::zeros); //hazard curves for infected to disease dead state
		
		cube state(5,days.n_elem,gamma1.n_elem,fill::zeros); //container to monitor individuals' states
			//row 1 = susceptible, row 2 = dead, row 3 = infected, row4 = recovered, row5 = disease death
			//columns - day
			//stack - MCMC rep
			
	
		//Create survival curves for each state - for each MCMC rep
		#pragma omp parallel for
		for(uword i=0; i<gamma1.n_elem; ++i){ 
				double UCH1 = 0; //Unit cumulative hazard
				double CH01 = 0; //Cumulative hazard
				double UCH2 = 0;
				double CH02 = 0;
				double UCH3 = 0;
				double CH03 = 0;
				double UCH4 = 0;
				double CH04 = 0;
				double UCH6 = 0;
				double CH06 = 0;
				double pos2 = 1/(1+exp(-beta01(i,0))); //prob infected @ 400 strain, @ mean elisa value, @ mean distance -last 2=0 since standardized covariates
				double p0i1 = exp(beta0dd(i,0))/(1 + exp(beta0dd(i,0)) + exp(beta0rec(i,0)) );//1/(1+exp(-beta0dd(i,0))); //prob disease death @ 400 strain, mean elisa value
				double p0i2 = exp(beta0rec(i,0))/(1 + exp(beta0dd(i,0)) + exp(beta0rec(i,0)) );//1/(1+exp(-beta0rec(i,0))); //prob recovered @ 400 strain, mean elisa value
				double p0r = 1/(1+exp(-beta06(i,0))); //prob reinfection 
			for(uword j=0; j<Ka.n_cols; ++j){
				
				
				UCH1 = exp(gamma1(i));
				CH01 = CH01 + UCH1;
				UCH2 = exp(gamma2(i) + Ka(i,j));
				CH02 = CH02 + UCH2;
				
				Ssucept(i,j) = pos2*exp(-CH02)+(1-pos2)*exp(-CH01);  
				hazinf(i,j) = UCH2; //hazard
			}
			    UCH1 = 0; //Unit cumulative hazard
				CH01 = 0; //Cumulative hazard
			
			for(uword j=0; j<Kadd.n_cols; ++j){
				UCH1 = exp(gamma1(i));
				CH01 = CH01 + UCH1;	
				
				UCH3 = exp(gamma3(i)+Kadd(i,j));
				CH03 = CH03 + UCH3;
				UCH4 = exp(gamma4(i));
				CH04 = CH04 + UCH4;
				UCH6 = exp(gamma6(i));
				CH06 = CH06 + UCH6;
			
				Sinfect(i,j) = p0i1*exp(-CH03)+p0i2*exp(-CH04)+(1-p0i1-p0i2)*exp(-CH01);
				Srec(i,j) = p0r*exp(-CH06)+(1-p0r)*exp(-CH01);
				
				hazdd(i,j) = UCH3; //hazard
				
			}
		}
		

		#pragma omp parallel for
		for(uword i=0; i<gamma1.n_elem; ++i){ 
				//Draw covariates
				mat X(N,XX.n_cols,fill::zeros);
				for(uword w=0; w<X.n_rows; ++w){
					double draw=as_scalar(round(randu(1)*(XX.n_rows-1)));
					X(w,span::all) = XX(draw,span::all);
				}
				
				uvec index6b = {0, 2, 3, 4, 5};  //indices for covariates for recovery probability
 				
				for(int n=0; n<N; ++n){
				
					uvec index6a(1);
					index6a(0) = as_scalar(n);
					
					uvec index2b = {0, 1, 6, 9, 2, 3, 4, 5}; //index for prob 2 - covariates
				
					//Rcpp::Rcout<<"BEGIN"<< n<<std::endl; 
						
					//Transition probabilities
					vec probs(5,fill::zeros);  //probability of going into each ABSORPTIVE state
					
					
					//Susceptible
					double iprob = as_scalar(1/(1+exp(-X(n,span(0,5))*trans(beta01(i,span::all))))); //prob suscept to infected
					probs(0) = 1-iprob; //prob suscept to dead
					
					//Infected
					probs(1) = as_scalar(iprob*(exp(X.submat(index6a,index2b)*trans(beta0dd(i,span::all)))/
						(1+exp(X.submat(index6a,index2b)*trans(beta0dd(i,span::all))) + exp(X.submat(index6a,index2b)*trans(beta0rec(i,span::all)))))); //prob infected to disease death
						
					double rprob = as_scalar(iprob*(exp(X.submat(index6a,index2b)*trans(beta0rec(i,span::all)))/
						(1+exp(X.submat(index6a,index2b)*trans(beta0dd(i,span::all))) + exp(X.submat(index6a,index2b)*trans(beta0rec(i,span::all)))))); //prob infected to recovered
					
					probs(2) = as_scalar(iprob*(1-(probs(1)/iprob + rprob/iprob))); //prob infected to dying
					
					//Recovered
					double riprob = as_scalar(iprob*rprob*(1/(1+exp(-X.submat(index6a,index6b)*trans(beta06(i,span::all)))))); //prob recovered to infected
					probs(3) = as_scalar(iprob*rprob*(1-riprob/(iprob*rprob))); //prob recovered to dead
					probs(4) = as_scalar(iprob*rprob*riprob*(probs(1)/iprob)); //reinfected to disease dead
					
					
					vec sdead(Ka.n_cols,fill::zeros);
					double CH0a = 0;
					
					for(uword j=0; j<Ka.n_cols; ++j){
						double UCHa = exp(gamma1(i));
						CH0a = CH0a + UCHa;
						sdead(j) = exp(-CH0a);				
					
					}
			
		
					int outcome = ProbSampleReplace(probs);  //final fate
					
					//Rcpp::Rcout<<mean(state,2)<<std::endl; 
					
					
					
					//Outcome 0 - S to D
					//Outcome 1 - S to I to DD
					//Outcome 2 - S to I to D
					//Outcome 3 - S to I to R to D
					//Outcome 4 - S to I to R to I to DD
					
					int dtime =0;
					int dtimea=0;
					int irtime =0;
					int itime =0;
					int rritime = 0;
					int ddtime = 0;
					
					double d1, d2, d3, d4, d5, d6;
					
					d1 = fail(sdead); //prob failing in each interval - non disease death
		
											
					//Susceptible Population
					//Determine timing of final fate and assign states for all days
					if(outcome == 0){ //die as susceptible		
						if(d1 >= days(days.n_elem-1)){
							state(span(0),span::all,span(i)) = state(span(0),span::all,span(i))+1; //survives entire time no infection
						}else{
							dtime = as_scalar(find(days>d1,1))-1; //determine which interval failed 
							state(span(0),span(0,dtime),span(i)) = state(span(0),span(0,dtime),span(i))+1; //set susceptible intervals
							state(span(1),span(dtime+1,days.n_elem-1),span(i)) = state(span(1),span(dtime+1,days.n_elem-1),span(i))+1; //set dead intervals
							}					
					 continue;	
					}
					
					//Infected population
					vec sinf(Ka.n_cols,fill::zeros);  //infection "survival"
					double CH0b = 0;
					uvec indices;
					indices << 1 << 7; 
					uvec indicesa(1);
					indicesa(0) = n;
					double eff1 = as_scalar(X.submat(indicesa,indices)*trans(betah(i,span::all))); //covariate effects 
					
					for(uword j=0; j<Ka.n_cols; ++j){
						double UCHb = exp(gamma2(i) + eff1 + Ka(i,j));
						CH0b = CH0b + UCHb;
						sinf(j) = exp(-CH0b);				
					}

					d2 = fail(sinf); //prob failing in each interval
					
					if(d2 >= days(days.n_elem-1)){ //infection after study end
							state(span(0),span::all,span(i)) = state(span(0),span::all,span(i))+1; //survives entire time no infection
							continue;
						}else{
							itime = as_scalar(find(days>d2,1))-1; //determine which interval failed 	 
							state(span(0),span(0,itime),span(i)) = state(span(0),span(0,itime),span(i))+1; //susceptible until infected
							state(span(2),span(itime+1,days.n_elem-1),span(i)) = state(span(2),span(itime+1,days.n_elem-1),span(i))+1; //infection span
						}
									
					
					if(outcome ==1){  //infected to disease death
					
						vec sdd(Kadd.n_cols,fill::zeros);  //infection "survival"
						double CH0c = 0;
						uvec indices2;
						indices2 << 1 << 6 ;
										
						uvec indices2a(1);
						indices2a(0) = n; 
						
						double eff2 = as_scalar(X.submat(indices2a,indices2)*trans(betah3(i,span::all))); //covariate effects 
						for(uword j=0; j<Kadd.n_cols; ++j){
							
							double UCHc = exp(gamma3(i) + eff2 + Kadd(i,j));
							CH0c = CH0c + UCHc;
							sdd(j) = exp(-CH0c);				
						}

						d3 = fail(sdd); //prob failing in each interval						
						double tot3 = d3+d2;
							
						
						if(tot3 >= days(days.n_elem-1)){ //if die of disease after study
								continue;
						}else{ //die of disease
							ddtime = as_scalar(find(days>tot3,1)); //determine which interval failed 	
//Rcpp::Rcout<<tot3<<":"<<n<<std::endl; 
							state(span(2),span(ddtime,days.n_elem-1),span(i)) = state(span(2),span(ddtime,days.n_elem-1),span(i))-1;
							state(span(4),span(ddtime,days.n_elem-1),span(i)) = state(span(4),span(ddtime,days.n_elem-1),span(i))+1;
							continue;	
							
						}	
					}
//Rcpp::Rcout<<":"<<n<<std::endl; 					
					if(outcome==2){//infected to non-disease death
						double tot4 = d1+d2;
						if(tot4 >= days(days.n_elem-1)){ //if die of disease after study
								continue;
						}else{ //die
							int idtime = as_scalar(find(days>tot4,1)); //determine which interval failed 
							state(span(2),span(idtime,days.n_elem-1),span(i))= state(span(2),span(idtime,days.n_elem-1),span(i))-1;
							state(span(1),span(idtime,days.n_elem-1),span(i))= state(span(1),span(idtime,days.n_elem-1),span(i))+1;
							continue;							
						}					   
					}
									
				
					if(outcome ==3 || outcome==4){ //infected to recovered
						vec srec(Kadd.n_cols,fill::zeros);  //infection "survival"
						double CH0d = 0;
						uvec indices3; 
						indices3 << 1 << 6 << 2 << 3 << 4 << 5; 
						uvec indices3a(1);
						indices3a(0) = n;
						double eff3 = as_scalar(X.submat(indices3a,indices3)*trans(betah4(i,span::all))); //covariate effects 
					
						for(uword j=0; j<Kadd.n_cols; ++j){
							double UCHd = exp(gamma4(i) + eff3);
							CH0d = CH0d + UCHd;
							srec(j) = exp(-CH0d);				
						}
						
						d4 = fail(srec); //prob failing in each interval
						double tot5 = d4+d2;				
								

						if(tot5 >= days(days.n_elem-1)){ //if recover after study
								continue;
						}else{ //recover
							irtime = as_scalar(find(days>tot5,1)); //determine which interval failed 
							state(span(2),span(irtime,days.n_elem-1),span(i))= state(span(2),span(irtime,days.n_elem-1),span(i))-1;
							state(span(3),span(irtime,days.n_elem-1),span(i))= state(span(3),span(irtime,days.n_elem-1),span(i))+1;
						}
					}
					
					if(outcome==3){ //recovered to non-disease death
						double tot6=d1+d2+d4;
						
						if(tot6 >= days(days.n_elem-1)){ //if die after study
								continue;
						}else{ //die 
							double rdtime = as_scalar(find(days>tot6,1)); //determine which interval failed 
							state(span(3),span(rdtime,days.n_elem-1),span(i))= state(span(3),span(rdtime,days.n_elem-1),span(i))-1;
							state(span(1),span(rdtime,days.n_elem-1),span(i))= state(span(1),span(rdtime,days.n_elem-1),span(i))+1;
							continue;	
						}			
							
					}
									
				
					if(outcome==4){ //recovered to re-infection
						vec sreinf(Kadd.n_cols,fill::zeros);  //infection "survival"
						double CH0e = 0;						
					
						for(uword j=0; j<Kadd.n_cols; ++j){
							double UCHe = exp(gamma6(i));
							CH0e = CH0e + UCHe;
							sreinf(j) = exp(-CH0e);				
						}
						
						d5 = fail(sreinf); //prob failing in each interval	
						double tot7	= d5+d4+d2;	
						
						
						if(tot7 >= days(days.n_elem-1)){ //if recovered to end of study
								continue;
						}else{ //reinfect
							rritime = as_scalar(find(days>tot7,1)); //determine which interval failed 
							state(span(3),span(rritime,days.n_elem-1),span(i))= state(span(3),span(rritime,days.n_elem-1),span(i))-1;
							state(span(2),span(rritime,days.n_elem-1),span(i))= state(span(2),span(rritime,days.n_elem-1),span(i))+1;
						}
						
						//reinfection to disease death
						vec sdd(Kadd.n_cols,fill::zeros);  //infection "survival"
						double CH0c = 0;
						uvec indices2;
						indices2 << 1 << 6 ;
						uvec indices2a(1);
						indices2a(0) = n; 
						double eff2 = as_scalar(X.submat(indices2a,indices2)*trans(betah3(i,span::all))); //covariate effects 
					
						for(uword j=0; j<Kadd.n_cols; ++j){
							double UCHc = exp(gamma3(i) + eff2 + Kadd(i,j));
							CH0c = CH0c + UCHc;
							sdd(j) = exp(-CH0c);				
						}
						
						d6 = fail(sdd); //prob failing in each interval);
						double tot8 = d2+d4+d5+d6;
						 
						
						if(tot8 >= days(days.n_elem-1)){ //if die of disease after study
								continue;
						}else{ //die of disease
							int ddtime2 = as_scalar(find(days>tot8,1)); //determine which interval failed
							state(span(2),span(ddtime2,days.n_elem-1),span(i))= state(span(2),span(ddtime2,days.n_elem-1),span(i))-1;
							state(span(4),span(ddtime2,days.n_elem-1),span(i))= state(span(4),span(ddtime2,days.n_elem-1),span(i))+1;
							continue;	
						}
			
					
					}

					
				}
					
		}

		
		//Summarize population values - infected
		mat muout = mean(state,2);
		mat temp025(state.n_rows,state.n_cols,fill::zeros);
		mat temp975(state.n_rows,state.n_cols,fill::zeros);		
		for(uword ii=0; ii<state.n_rows; ++ii){
			for(uword jj=0; jj<state.n_cols; ++jj){
				temp025(ii,jj) = quantilecpp(state(span(ii),span(jj),span::all),0.025);
				temp975(ii,jj) = quantilecpp(state(span(ii),span(jj),span::all),0.975);
				
			}
		}

		//Summarize survival curves
		mat Ss(Ssucept.n_cols,3,fill::zeros);
		mat Si(Sinfect.n_cols,3,fill::zeros); 
		mat Sr(Srec.n_cols,3,fill::zeros);  
		mat hi(hazinf.n_cols,3,fill::zeros);
		mat hdd(hazdd.n_cols,3,fill::zeros);
		
		Ss(span::all,0) = trans(median(Ssucept,0));
		Si(span::all,0) = trans(median(Sinfect,0));
		Sr(span::all,0) = trans(median(Srec,0));
		
		hi(span::all,0) = trans(median(hazinf,0));
		hdd(span::all,0) = trans(median(hazdd,0));
		
		#pragma omp parallel for
		for(uword kk=0; kk<Ssucept.n_cols; ++kk){
			Ss(kk,1) = quantilecpp(Ssucept(span::all,kk),0.025);
			Ss(kk,2) = quantilecpp(Ssucept(span::all,kk),0.975);
		}
		
		#pragma omp parallel for
		for(uword kk=0; kk<Sinfect.n_cols; ++kk){
			Si(kk,1) = quantilecpp(Sinfect(span::all,kk),0.025);
			Si(kk,2) = quantilecpp(Sinfect(span::all,kk),0.975);
		}
		
		#pragma omp parallel for	
		for(uword kk=0; kk<Srec.n_cols; ++kk){
			Sr(kk,1) = quantilecpp(Srec(span::all,kk),0.025);
			Sr(kk,2) = quantilecpp(Srec(span::all,kk),0.975);
		}
		
		#pragma omp parallel for
		for(uword kk=0; kk<hazinf.n_cols; ++kk){
			hi(kk,1) = quantilecpp(hazinf(span::all,kk),0.025);
			hi(kk,2) = quantilecpp(hazinf(span::all,kk),0.975);
		}
		
		#pragma omp parallel for
		for(uword kk=0; kk<hazdd.n_cols; ++kk){
			hdd(kk,1) = quantilecpp(hazdd(span::all,kk),0.025);
			hdd(kk,2) = quantilecpp(hazdd(span::all,kk),0.975);
		}			
				
		return List::create(Rcpp::Named("Nmu") = muout,
							Rcpp::Named("quantile025") = temp025,
							Rcpp::Named("quantile975") = temp975,
							Rcpp::Named("Ssucept") = Ss,
							Rcpp::Named("Sinfect") = Si,
							Rcpp::Named("Srec") = Sr,
							Rcpp::Named("hazinf") = hi,
							Rcpp::Named("hazdd") = hdd
							);
}			


 

