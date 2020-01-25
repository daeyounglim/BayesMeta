#ifndef BAYESMETA_DICNMR_H
#define BAYESMETA_DICNMR_H
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo, RcppNumerical)]]
using namespace Numer;
using namespace arma;

class ObservedLik: public Func
{
public:
	ObservedLik(int& Tk_,
				arma::vec& resid_,
			    double& lam_k_,
			    arma::vec& Z_k_,
			    arma::mat& E_k_,
			    arma::vec& sig2_k_,
			    arma::vec& npt_k_,
			    arma::mat& Rho_,
			    double& nu_) : Tk(Tk_), resid(resid_), lam_k(lam_k_), Z_k_(Z_k_), E_k(E_k_), sig2_k_(sig2_k), npt_k(npt_k_), Rho(Rho_), nu(nu) {}

	double operator()(const double& lam) const {
		double loglik = -M_LN_SQRT_2PI * static_cast<double>(Tk) + 
			0.5 * nu * (std::log(nu) - M_LN2) - R::lgammafn(0.5 * nu) + (0.5 * nu - 1.0) * lam - 0.5 * nu * lam;
		double logdet_val;
		double logdet_sign;
		mat ERE_S = diagmat(Z_k) * E_k.t() * Rho * E_k * diagmat(Z_k / lam);
		ERE_S.diag() += sig2_k / npt_k;
		log_det(logdet_val, logdet_sign, ERE_S);
		return std::exp(loglik -0.5 * (logdet_val + arma::accu(resid % arma::solve(ERE_S, resid))));
	}
};

double int_ObservedLik(int& Tk,
					   arma::vec& resid,
			    	   double& lam_k,
			    	   arma::vec& Z_k,
			    	   arma::mat& E_k,
			    	   arma::vec& sig2_k,
			    	   arma::vec& npt_k,
			    	   arma::mat& Rho,
			    	   double& nu);

#endif
