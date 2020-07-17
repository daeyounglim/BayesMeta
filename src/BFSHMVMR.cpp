#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "ListBuilder.h"
#include "misc_nmr.h"
#include "nelmin.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

/*
C++ code for fitting
Bayesian Flexible Hierarchical Skew Heavy-Tailed Multivariate
Meta Regression Models for Individual Patient Data

x    : `xcols`-dimensional
beta : xcols * J
gamma: q * J
delta_star : J * K

*/


// [[Rcpp::export]]
Rcpp::List BFSHMVMR(const arma::mat& Outcome,
					const arma::mat& Covariate,
					const arma::vec& Treat,
					const arma::uvec& Patient,
					const arma::uvec& Trial,
					const arma::uvec& Npt,
					const double& a1,
					const double& a2,
					const double& a3,
					const double& a4,
					const double& a5,
					const double& b1,
					const double& b2,
					const double& b3,
					const double& b4,
					const double& b5,
					const double& c1,
					const double& c2,
					const double& c3,
					const double& d1,
					const double& d2,
					const int& N, // number of 
					const int& ndiscard,
					const int& nskip,
					const int& nkeep,
					const bool verbose
					) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	const int N = Outcome.n_rows;
	const int J = Outcome.n_cols;
	const int K = Npt.n_elem;
	const int xcols = Covariate.n_cols;
	const int nw = 2;
	const int nt = (xcols + nw) * J;
	const int sum_nk = arma::accu(Npt);


	arma::field<arma::uvec> idxks(K);
	for (int k = 0; k < K; ++k) {
		uvec idx = find(ids == k+1);
		idxks(k) = idx;
	}

	vec theta(nt, fill::zeros);
	vec lambda(N, fill::ones);
	vec delta(J, fill::zeros);
	vec phi(N, fill::ones); // `zeta` in FORTRAN code
	mat zz(N, J, fill::ones);
	vec Ez(J, fill::ones);
	mat Sigma(J, J, fill::eye);
	mat Sigmainv = Sigma;
	cube Omega(2,2,J, fill::zeros);
	for (int j = 0; j < J; ++j) {
		Omega.slice(j) = eye<mat>(2, 2);
	}
	mat gamR(nw*J, K, fill::zeros); // `xi` in FORTRAN code

	double nu = 1.0; // `vv` in FORTRAN code
	double nu0 = 1.0; // `v0` in FORTRAN code
	double sig2_theta = 1.0;
	double sig2_delta = 1.0;

	/*******************
	Begin burn-in period
	*******************/
	if (verbose) {
		Rcout << "Warming up" << endl;
	}
	{
		Progress prog(ndiscard, verbose);
		for (int idiscard = 0; idiscard < ndiscard; ++idiscard) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
			/********
			Update zz
			********/
			mat DS = arma::diagmat(delta) * Sigmainv;
			mat DSD = DS * arma::diagmat(delta);
			for (int i = 0; i < N; ++i) {
				rowvec x_i = Covariate.row(i);
				rowvec w_i(2, fill::ones);
				w_i(1) = Treat(i);
				mat X(J, xcols * J);
				mat W(J, 2*J);
				for (int j = 0; j < J; ++j) {
					X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
					W(j, span(j*2, (j+1)*2-1)) = w_i;
				}
				mat Az = lambda(i) * DSD;
				mat Bz = lambda(i) * DS * 
			}
		}
	}

	return ListBuilder()
	.add("beta", 0.0);
}