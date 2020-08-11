#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "random.h"
#include "linearalgebra.h"
#include "ListBuilder.h"
#include "nelmin.h"

// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

/*
// #include "misc_nmr.h"
C++ code for fitting
Bayesian multivariate skew meta-regression models for individual patient data

x    : `xcols`-dimensional
beta : xcols * J
gamma: q * J
delta_star : J * K


If J = 1, don't sample psi
*/

// [[Rcpp::export]]
Rcpp::List BMVMR_POCov(const arma::mat& Outcome,
					   const arma::mat& SD,
					   const arma::mat& XCovariate,
					   const arma::mat& WCovariate,
					   const arma::vec& Treat,
					   const arma::uvec& Trial,
					   const arma::vec& Npt,
					   const double& c0,
					   const double& dj0, // hyperparameter for Omega
					   const double& d0, // hyperparameter for Sigma
					   const double& s0,
					   const arma::mat& Omega0,
					   const arma::mat& Sigma0,
					   const int& K, // # of Trials
					   const int& T, // # of Treatments
					   const int& fmodel,
					   const int& ndiscard,
					   const int& nskip,
					   const int& nkeep,
					   const bool& verbose) {
	using namespace arma;
	using namespace std;
	using namespace Rcpp;
	using namespace R;

	const int N = Outcome.n_rows;
	const int J = Outcome.n_cols;
	const int xcols = XCovariate.n_cols;
	const int nw = WCovariate.n_cols;
	const int nt = (xcols + nw) * J;

	arma::field<arma::uvec> idxks(K);
	// mat n_tk(T, K, fill::zeros);
	for (int k = 0; k < K; ++k) {
		uvec idx = find(Trial == k);
		idxks(k) = idx;
	}


	/***********************
	Parameter Initialization
	***********************/

	vec theta(nt, fill::zeros);
	mat gamR(nw*J, K, fill::ones);
	mat Omegainv(nw*J, nw*J, fill::eye);
	arma::field<arma::mat> Sigmainvs(T,K);
	mat Rtk(N, J * (J - 1) / 2, fill::zeros);
	for (int k = 0; k < K; ++k) {
		for (int t = 0; t < T; ++t) {
			Sigmainvs(t,k) = eye<mat>(J,J);
		}
	}
	mat Rho(J, J, fill::eye);

	/************
	Miscellaneous
	************/
	const mat Omega0inv = arma::inv(Omega0);
	const mat Sigma0inv = arma::inv(Sigma0);
	const df = s0 + arma::accu(Npt);
	const shape_omega = static_cast<double>(K) + dj0;
	mat resid = Outcome;
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


			/***********
			Update theta
			***********/
			mat Sig_theta(nt, nt, fill::zeros);
			Sig_theta.diag().fill(1.0 / c01);
			vec mu_theta(nt, fill::zeros);
			for (int k = 0; k < K; ++k) {
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				mat XSX(ntd, ntd, fill::zeros);
				mat WSX(nw*J, ntd, fill::zeros);
				vec WSy(nw*J, fill::zeros);
				vec XSy(ntd, fill::zeros);
				mat Sig_gamk = Omegainv;
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = XCovariate.row(i_k);
					rowvec w_i = WCovariate.row(i_k);
					rowvec y_i = Outcome.row(i_k);
					int t = Treat(i_k);
					mat Sigmainv = Sigmainvs(t, k);
					double ntk = Npt(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					XSX += ntk * Xstar.t() * Sigmainv * Xstar;
					XSy += ntk * Xstar.t() * Sigmainv * y_i.t();
					mat WS = ntk * W.t() * Sigmainv;
					Sig_gamk += WS * W;
					WSX += WS * Xstar;
					WSy += WS * y_i.t();
				}
				mat Sig_gamk_inv = arma::inv(Sig_gamk);
				Sig_theta += XSX - WSX.t() * Sig_gamk_inv * WSX;
				mu_theta += XSy - WSX.t() * Sig_gamk_inv * WSy;
			}
			Sig_theta = 0.5 * (Sig_theta + Sig_theta.t());
			mat Sig_thetaChol = chol(Sig_theta);
			vec atmp(nt);
			std::generate(atmp.begin(), atmp.end(), ::norm_rand);
			theta = arma::solve(arma::trimatu(Sig_thetaChol), arma::solve(trimatl(Sig_thetaChol.t()), mu_theta) + atmp);
			for (int k = 0; k < K; ++k) {
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = XCovariate.row(i_k);
					rowvec w_i = WCovariate.row(i_k);
					rowvec y_i = Outcome.row(i_k);
					int t = Treat(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					vec resid_i = y_i.t() - join_horiz(X, W) * theta;
					resid.row(i_k) = resid_i.t();
				}
			}


			/**********
			Update gamR
			**********/
			for (int k = 0; k < K; ++k) {
				mat Siggam = Omegainv;
				vec mugam(nw*J, fill::zeros);
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = XCovariate.row(i_k);
					rowvec w_i = WCovariate.row(i_k);
					rowvec y_i = Outcome.row(i_k);
					double ntk = Npt(i_k);
					int t = Treat(i_k);
					mat Sigmainv = Sigmainvs(t, k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					mat WS = W.t() * Sigmainv;
					Siggam += ntk * WS * W;
					mugam += ntk * WS * (y_i.t() - arma::join_horiz(X,W) * theta);
				}
				Siggam = 0.5 * (Siggam + Siggam.t());
				mat SiggamChol = chol(Siggam);
				vec gtmp(nw*J);
				std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
				gamR.col(k) = arma::solve(arma::trimatu(SiggamChol), arma::solve(trimatl(SiggamChol.t()), mugam) + gtmp);
			}

			/***********
			Update Omega
			***********/
			for (int jj = 0; jj < J; ++jj) {
				mat qq = S1inv;
				vec gamstar(nw);
				for (int k = 0; k < K; ++k) {
					for (int j = 0; j < nw; ++j) {
						gamstar(j) = gamR(nw*jj+j, k);
					}

					for (int j1 = 0; j1 < nw; ++j1) {
						for (int j2 = 0; j2 < nw; ++j2) {
							qq(j1, j2) += gamstar(j1) * gamstar(j2);
						}
					}
				}
				Omegainv(span(jj*nw, (jj+1)*nw-1), span(jj*nw, (jj+1)*nw-1)) = rwish(shape_omega, qq.i());
			}

			/***********
			Update Sigma
			***********/
			if (fmodel == 1) {
				for (int i = 0; i < N; ++i) {
					int k = Trial(i);
					int t = Treat(i);
					double ntk = Npt(i);
					double shape = d0 + 0.5 * ntk;
					rowvec sd2_i = arma::square(SD.row(i));
					vec gam_k = gamR.col(k);
					rowvec w_i = WCovariate.row(i);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					vec resid_i2 = arma::square(arma::trans(resid.row(i)) - W * gam_k);
					vec sig2inv(J);
					for (int j = 0; j < J; ++j) {
						double rate = s0 + ntk * 0.5 * resid_i2(j) + (ntk - 1.0) * 0.5 * sd2_i(j);
						sig2inv(j) = ::Rf_rgamma(shape, 1.0) / rate;
					}
					Sigmainvs(t, k) = arma::diagmat(sig2inv);
				}
			} else {
				if (fmodel == 2) {
					mat qq = Sigma0inv;
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						int t = Treat(i);
						double ntk = Npt(i);
						vec gam_k = gamR.col(k);
						mat R = veclinv(trans(Rtk.row(i)));
						R.diag().fill(1.0);
						rowvec w_i = WCovariate.row(i);
						mat V = arma::diagmat(SD.row(i));
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
						qq += ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;
					}
					mat Siginv_new = rwish(df, qq);
				} else if (fmodel == 3) {

				}

				// Update R
				
			}
		}

	}
}