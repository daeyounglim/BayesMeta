#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <Rdefines.h>
#include "neldermead.h"
#include "dic_nmr.h"
// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppNumerical)]]
using namespace Numer;
using namespace arma;


/********************************
negative log-likelihood of lambda
for Nelder-Mead method
********************************/
double neg_loglik_lam(const double& lam_k,
				  	  const double& nu,
				  	  const arma::vec& resid_k, // = y_k - X_k * beta
				  	  const arma::vec& Z_k,
				  	  const arma::mat& ERE, // = E_k' * Rho * E_k
				  	  const arma::vec& sig2_k) {
	using namespace arma;

	double eta_k = std::log(lam_k);
	double loglik = 0.5 * nu * (eta_k - lam_k) - eta_k;
	mat tmpmat = arma::diagmat(Z_k) * ERE * arma::diagmat(Z_k / lam_k);
	tmpmat.diag() += sig2_k;
	double logdet_val;
	double logdet_sign;
	log_det(logdet_val, logdet_sign, tmpmat);
	loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(tmpmat, resid_k));
	return -loglik;
})

double int_ObservedLik(int& Tk,
					   arma::vec& resid,
			    	   double& lam_k,
			    	   arma::vec& Z_k,
			    	   arma::mat& E_k,
			    	   arma::vec& sig2_k,
			    	   arma::vec& npt_k,
			    	   arma::mat& Rho,
			    	   double& nu,
			    	   double& maxll) {
	ObservedLik f(Tk, resid, lam_k, Z_k, E_k, sig2_k, npt_k, Rho, nu, maxll);
	double err_est;
	int err_code;
	return Numer::integrate(f, 0.0, R_PosInf, err_est, err_code);
}

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
	vec phi_est = arma::mean(phis, 1);
	mat Rho_est = arma::mean(Rhos, 2);

	vec DICs(nkeep, fill::zeros);
	{
		Progress prog(nkeep, true);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}

			vec beta_ikeep = betas.col(ikeep);
			vec sig2_ikeep = sig2s.col(ikeep);
			vec phi_ikeep = phis.col(ikeep);
			mat Rho_ikeep = Rhos.slice(ikeep);
			vec Z_ikeep = arma::exp(z * phi_ikeep);

			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				vec y_k = y(idx);
				mat E_k = Eks(k);
				mat X_k = arma::join_horiz(Xks(k), E_k.t());
				vec resid_k = y_k - X_k * beta_ikeep;
				vec Z_k = Z_ikeep(idx);
				mat ERE = E_k.t() * Rho_ikeep * E_k;
				vec sig2_k = sig2_ikeep(idx) / npt(idx);

				vec start(1);
				start(0) = std::log(sig2_est(k));
				vec xmin(1, fill::zeros);
				double maxll = 0.0; // = ynewlo
				double reqmin = 1.0e-20;
				int konvge = 5;
				int kcount = 1000;
				vec step(1);
				step(0) = 0.2;
				int icount = 0;
				int numres = 0;
				int ifault = 0;

				nelmin ( double fn ( const arma::vec& x ),
			             1,
			             arma::vec& start,
			             arma::vec& xmin, 
			             double& ynewlo,
			             double& reqmin,
			             arma::vec& step,
			             int konvge,
			             int kcount, 
			             int& icount,
			             int& numres,
			             int& ifault )
				neg_loglik_lam(const double& lam_k,
				  	  const double& nu,
				  	  const arma::vec& resid_k, // = y_k - X_k * beta
				  	  const arma::vec& Z_k,
				  	  const arma::mat& ERE, // = E_k' * Rho * E_k
				  	  const arma::vec& sig2_k
			}

			prog.increment();
		}
	}
}















