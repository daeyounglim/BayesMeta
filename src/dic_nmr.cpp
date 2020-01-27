#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <RcppNumerical.h>
#include <Rdefines.h>
#include "dic_nmr.h"
// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppNumerical)]]
using namespace Numer;
using namespace arma;


double int_observedlik(const int& Tk,
					   const arma::vec& resid,
			    	   const arma::mat& ZEREZ_S,
			    	   const double& nu,
			    	   const double& maxll) {
	ObservedLik f(Tk, resid, ZEREZ_S, nu, maxll);
	double err_est;
	int err_code;
	return maxll + std::log(Numer::integrate(f, 0.0, R_PosInf, err_est, err_code, 100, 1.0e-08, 1.0e-06, Integrator<double>::GaussKronrod21));
}

/**************************************************
Calculate the goodness of fit measures

+ Dev(theta) = -2 * log L(theta | D_oy)
+ p_D = E[Dev(theta)] - Dev(thetabar)
+ DIC = Dev(thetabar) + 2 * p_D
***************************************************/

// [[Rcpp::export]]
Rcpp::List calc_modelfit(const arma::vec& y,
						 const arma::mat& x,
						 const arma::mat& z,
						 const arma::uvec& ids,
						 const arma::uvec& iarm,
						 const arma::vec& npt,
						 const double& nu,
						 const arma::mat& betas,
						 const arma::mat& sig2s,
						 const arma::mat& phis,
						 const arma::mat& lams,
						 const arma::cube& Rhos,
						 const int& K,
						 const int& nT,
						 const int& nkeep) {
	using namespace arma;

	/* make a list of y_k, X_k, z_k*/
	arma::field<arma::mat> Xks(K);
	arma::field<arma::mat> Eks(K);
	arma::field<arma::uvec> idxks(K);
	for (int k = 0; k < K; ++k) {
		uvec idx = find(ids == k+1);
		idxks(k) = idx;
		Xks(k) = x.rows(idx);
		int idx_l = idx.n_elem;
		mat Ek(nT, idx_l, fill::zeros);
		
		uvec iarm_k = iarm(idx);
		for (int j = 0; j < idx_l; ++j) {
			Ek(iarm_k(j),j) = 1.0;
		}
		Eks(k) = Ek;
	}

	vec beta_est = arma::mean(betas, 1);
	vec sig2_est = arma::mean(sig2s, 1);
	vec lam_est = arma::mean(lams, 1);
	vec phi_est = arma::mean(phis, 1);
	mat Rho_est(nT, nT, fill::zeros);
	for (int ikeep=0; ikeep < nkeep; ++ikeep) {
		Rho_est += Rhos.slice(ikeep);
	}
	Rho_est /= static_cast<double>(nkeep);
	// mat Rho_est = arma::mean(Rhos, 2);

	/*******************************************
	Dev(thetabar) = -2 * log L(thetabar | D_oy)
	*******************************************/
	double Dev_thetabar = 0.0;
	vec Z_est = arma::exp(z * phi_est);
	for (int k=0; k < K; ++k) {
		uvec idx = idxks(k);
		vec y_k = y(idx);
		mat E_k = Eks(k);
		mat X_k = arma::join_horiz(Xks(k), E_k.t());
		vec resid_k = y_k - X_k * beta_est;
		vec Z_k = Z_est(idx);
		mat ZEREZ_S = diagmat(Z_k) * E_k.t() * Rho_est * E_k * diagmat(Z_k / lam_est(k));
		vec sig2_k = sig2_est(idx) / npt(idx);
		ZEREZ_S.diag() += sig2_k;

		NegLogLikLam nlll(nu, resid_k, ZEREZ_S);
		nlll.optimize(std::log(lam_est(k)));
		double maxll = -nlll.ynewlo;

		int Tk = idx.n_elem;
		Dev_thetabar += -2.0 * int_observedlik(Tk, resid_k, ZEREZ_S, nu, maxll);
	}



	double Dev_bar = 0.0;
	{
		Progress prog(nkeep, true);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}

			vec beta_ikeep = betas.col(ikeep);
			vec sig2_ikeep = sig2s.col(ikeep);
			vec phi_ikeep = phis.col(ikeep);
			vec lam_ikeep = lams.col(ikeep);
			mat Rho_ikeep = Rhos.slice(ikeep);
			vec Z_ikeep = arma::exp(z * phi_ikeep);

			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				vec y_k = y(idx);
				mat E_k = Eks(k);
				mat X_k = arma::join_horiz(Xks(k), E_k.t());
				vec resid_k = y_k - X_k * beta_ikeep;
				vec Z_k = Z_ikeep(idx);
				double lam_k = lam_ikeep(k);
				mat ZEREZ_S = diagmat(Z_k) * E_k.t() * Rho_ikeep * E_k * diagmat(Z_k / lam_k);
				vec sig2_k = sig2_ikeep(idx) / npt(idx);
				ZEREZ_S.diag() += sig2_k;

				NegLogLikLam nlll(nu, resid_k, ZEREZ_S);
				nlll.optimize(std::log(lam_k));
				double maxll = -nlll.ynewlo;

				int Tk = idx.n_elem;
				Dev_bar += -2.0 * int_observedlik(Tk, resid_k, ZEREZ_S, nu, maxll);
			}

			prog.increment();
		}
		Dev_bar /= static_cast<double>(nkeep);
		double p_D = Dev_bar - Dev_thetabar;
		double DIC = Dev_thetabar + 2.0 * p_D;
		return Rcpp::List::create(Rcpp::Named("DIC")=DIC);
	}
}















