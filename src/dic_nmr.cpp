#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <Rdefines.h>
#include "dic_nmr.h"
// [[Rcpp::depends(RcppArmadillo, RcppNumerical)]]
using namespace Numer;
using namespace arma;

double int_ObservedLik(int& Tk,
					   arma::vec& resid,
			    	   double& lam_k,
			    	   arma::vec& Z_k,
			    	   arma::mat& E_k,
			    	   arma::vec& sig2_k,
			    	   arma::vec& npt_k,
			    	   arma::mat& Rho,
			    	   double& nu) {
	ObservedLik f(Tk, resid, lam_k, Z_k, E_k, sig2_k, npt_k, Rho, nu);
	double err_est;
	int err_code;
	return integrate(f, 0.0, R_PosInf, err_est, err_code);
}
