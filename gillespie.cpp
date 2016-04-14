// **************************************************
// WITS MODEL
// C++ CODE TO RUN Gillespie algorithm
// FOR USE WITH RCPP
// **************************************************

// Olivier Restif, version 1, May 2014

// This file contains all the C++ code to be exported to R using Rcpp
// Use of multiple source files require development of a package.

// ---------------- Notes ------------------------
// Random Numbers: 
//	- For non-prallelised code, simply use functions provided by Rcpp
// 	- For paralelised code, use the trng library which is compatible with OpenMP


#include <Rcpp.h>
using namespace Rcpp;


// Run one simulation from t_start to t_end with constant parameters and a choice of initial conditions
// par must be a vector with 8 doubles: cL, cS, eL, eS, kL, kS, rL, rS
// n_init must be a vector of integers nB,nL,nC
// Return vector of three integers: nB_out, nL_out, nS_out
// Use Rcpp RNG with seed set by R prior to call -- NOT SUITABLE FOR PARALLEL EXECUTION !

// [[Rcpp::export("single_wits_gillespie_r")]]
IntegerVector single_wits_gillespie_cpp(double t_start, double t_end, IntegerVector n_init, NumericVector par)
{
	double t=t_start, cL=par[0], cS=par[1], eL=par[2], eS=par[3], kL=par[4], kS=par[5], rL=par[6], rS=par[7], u=0.0;
	double rates[8];
	int nB=n_init[0], nL=n_init[1], nS=n_init[2], i=0;
	IntegerVector n_out(3);
	
	// Initialise Rcpp RNG
	RNGScope scope;
	
	// Start the loop
	while(t<t_end){
		// 1. Calculate cumulative rates
		rates[0] = cL*nB;
		rates[1] = cS*nB+rates[0];
		rates[2] = eL*nL+rates[1];
		rates[3] = eS*nS+rates[2];
		rates[4] = kL*nL+rates[3];
		rates[5] = kS*nS+rates[4];
		rates[6] = rL*nL+rates[5];
		rates[7] = rS*nS+rates[6];
		
		// Check for extinction: end function and return t
		if(rates[7]<1E-6){
			n_out[0] = nB;
			n_out[1] = nL;
			n_out[2] = nS;
			return n_out;
		}
		
		// 2. Draw the time to next event
		t += R::rexp(1/rates[7]);
		// Stop here if t > t_end
		if(t>t_end){
			n_out[0] = nB;
			n_out[1] = nL;
			n_out[2] = nS;
			return n_out;
		}
		// 3. Draw the next event
		u = R::runif(0,rates[7]);
		i=0;
		while(u>rates[i]) i++;
		
		// 4. Apply event i
		switch(i){
			case 0: nB--; nL++; break;
			case 1: nB--; nS++; break;
			case 2: nB++; nL--; break;
			case 3: nB++; nS--; break;
			case 4: nL--; break;
			case 5: nS--; break;
			case 6: nL++; break;
			case 7: nS++; break;
//			default: printf("\nError in choice of event.\n");
		}
		
	}
	n_out[0] = nB;
	n_out[1] = nL;
	n_out[2] = nS;
	return n_out;
}

// [[Rcpp::export("single_wits_gillespie_r_2")]]
NumericVector single_wits_gillespie_cpp_2(double t_start, double t_end, IntegerVector n_init, NumericVector par)
{
	double t=t_start, t_last=0.0, cL=par[0], cS=par[1], eL=par[2], eS=par[3], kL=par[4], kS=par[5], rL=par[6], rS=par[7], u=0.0;
	double rates[8];
	int nB=n_init[0], nL=n_init[1], nS=n_init[2], i=0;
	NumericVector n_out(4);
	
	// Initialise Rcpp RNG
	RNGScope scope;
	
	// Start the loop
	while(t<t_end){
		// 1. Calculate cumulative rates
		rates[0] = cL*nB;
		rates[1] = cS*nB+rates[0];
		rates[2] = eL*nL+rates[1];
		rates[3] = eS*nS+rates[2];
		rates[4] = kL*nL+rates[3];
		rates[5] = kS*nS+rates[4];
		rates[6] = rL*nL+rates[5];
		rates[7] = rS*nS+rates[6];
		
		// Check for extinction: end function and return 
		if(rates[7]<1E-6){
			n_out[0] = nB;
			n_out[1] = nL;
			n_out[2] = nS;
			n_out[3] = t;
			return n_out;
		}
		
		// 2. Draw the time to next event
		t_last = t;
		t += R::rexp(1/rates[7]);
		// Stop here if t > t_end
		if(t>t_end){
			n_out[0] = nB;
			n_out[1] = nL;
			n_out[2] = nS;
			n_out[3] = t_last;
			return n_out;
		}
		// 3. Draw the next event
		u = R::runif(0,rates[7]);
		i=0;
		while(u>rates[i]) i++;
		
		// 4. Apply event i
		switch(i){
		case 0: nB--; nL++; break;
		case 1: nB--; nS++; break;
		case 2: nB++; nL--; break;
		case 3: nB++; nS--; break;
		case 4: nL--; break;
		case 5: nS--; break;
		case 6: nL++; break;
		case 7: nS++; break;
		}
		
	}
	n_out[0] = nB;
	n_out[1] = nL;
	n_out[2] = nS;
	n_out[3] = t;
	return n_out;
}

