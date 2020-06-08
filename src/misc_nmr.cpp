#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rdefines.h>
#include "misc_nmr.h"
// [[Rcpp::depends(RcppArmadillo)]]

double loglik_lam(const double& eta_k,
				  const double& nu,
				  const arma::vec& resid_k, // = y_k - X_k * beta
				  const arma::mat& ERE_k, // = Z(phi) *  E_k' * Rho * E_k * Z(phi)
				  const arma::vec& sig2_k,
				  const int& Tk) { 
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	double lam_k = std::exp(eta_k);

	double loglik = (0.5 * nu - 1.0) * eta_k - 0.5 * nu * lam_k + 0.5 * nu * (std::log(nu) - M_LN2) - R::lgammafn(0.5 * nu);
	mat tmpmat = std::exp(-eta_k) * ERE_k;
	tmpmat.diag() += sig2_k;
	double logdet_val;
	double logdet_sign;
	log_det(logdet_val, logdet_sign, tmpmat);
	loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(tmpmat, resid_k)) - M_LN_SQRT_2PI * static_cast<double>(Tk);
	return loglik;
}

double loglik_eta(const double& eta_k,
				  const double& nu,
				  const arma::vec& resid_k, // = y_k - X_k * beta
				  const arma::vec& Z_k,
				  const arma::mat& ERE, // = E_k' * Rho * E_k
				  const arma::vec& sig2_k) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	double loglik = 0.5 * nu * (eta_k - std::exp(eta_k));
	mat tmpmat = std::exp(-eta_k) * arma::diagmat(Z_k) * ERE * arma::diagmat(Z_k);
	tmpmat.diag() += sig2_k;
	double logdet_val;
	double logdet_sign;
	log_det(logdet_val, logdet_sign, tmpmat);
	loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(tmpmat, resid_k));
	return loglik;
}

double loglik_phi(const arma::vec& phi,
				  const arma::mat& z,
				  const double& c02inv,
				  const arma::vec& lam,
				  const arma::vec& sig2,
				  const arma::mat& Rho,
				  const arma::vec& resid,
				  const arma::field<arma::mat>& Eks,
				  const arma::field<arma::uvec>& idxks) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	vec Z = arma::exp(z * phi);
	const int K = Eks.n_elem;
	double loglik = -0.5 * c02inv * arma::accu(phi % phi);
	for (int k=0; k < K; ++k) {
		uvec idx_k = idxks(k);
		mat Z_k = arma::diagmat(Z(idx_k));
		vec sig2_k = sig2(idx_k);
		mat E_k = Eks(k);
		mat ERE = E_k.t() * Rho * E_k;
		vec resid_k = resid(idx_k);

		mat tmpmat = Z_k * ERE * Z_k / lam(k);
		tmpmat.diag() += sig2_k;
		double logdet_val;
		double logdet_sign;
		log_det(logdet_val, logdet_sign, tmpmat);

		loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(tmpmat, resid_k));
	}
	return loglik;
}


arma::mat pRho_to_Rho(arma::mat& pRho) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	const int nT = pRho.n_rows;
	mat Rho = pRho;

	for (int iL=2; iL < nT; ++iL) {
		mat rmat(iL-1,iL-1, fill::zeros);

		for (int iR=1; iR < nT-iL+1; ++iR) {
			vec rvec1(nT-2, fill::zeros);
			vec rvec3(nT-2, fill::zeros);

			for (int i1=1; i1 < iL; ++i1) {
				rvec1(i1-1) = Rho(iR-1, iR+i1-1);
				rvec3(i1-1) = Rho(iR+i1-1, iR+iL-1);
			}
			double rr11 = 0.0;
			double rr13 = 0.0;
			double rr33 = 0.0;

			for (int i1=1; i1 < iL; ++i1) {
				for (int j1=1; j1 < iL; ++j1) {
					rmat(i1-1,j1-1) = Rho(iR+i1-1, iR+j1-1);
				}
			}
			mat rmatinv = rmat.i();
			for (int i1=1; i1 < iL; ++i1) {
				for (int j1=1; j1 < iL; ++j1) {
					rr11 += rvec1(i1-1) * rmatinv(i1-1,j1-1) * rvec1(j1-1);
					rr13 += rvec1(i1-1) * rmatinv(i1-1,j1-1) * rvec3(j1-1);
					rr33 += rvec3(i1-1) * rmatinv(i1-1,j1-1) * rvec3(j1-1);
				}
			}
			Rho(iR-1, iR+iL-1) = rr13 + pRho(iR-1,iR+iL-1) * std::sqrt((1.0 - rr11) * (1.0 - rr33));
			Rho(iR+iL-1,iR-1) = Rho(iR-1, iR+iL-1);
		}
	}
	return Rho;
}

double loglik_z(const double& zprho,
				const int& index1,
				const int& index2,
				const arma::mat& pRho,
				const arma::vec& lam,
				const arma::vec& sig2,
				const arma::vec& Z,
				const arma::vec& resid,
				const arma::field<arma::mat>& Eks,
				const arma::field<arma::uvec>& idxks) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	mat tmppRho = pRho;
	const int nT = pRho.n_rows;
	const int K = Eks.n_elem;
	tmppRho(index1, index2) = (std::exp(2.0 * zprho) - 1.0) / (std::exp(2.0 * zprho) + 1.0);
	tmppRho(index2, index1) = tmppRho(index1, index2);

	mat tmpRho = pRho_to_Rho(tmppRho);

	double loglik = 0.0;
	for (int k=0; k < K; ++k) {
		uvec idx_k = idxks(k);
		mat Z_k = arma::diagmat(Z(idx_k));
		vec sig2_k = sig2(idx_k);
		mat E_k = Eks(k);
		mat ERE = E_k.t() * tmpRho * E_k;
		vec resid_k = resid(idx_k);

		mat tmpmat = Z_k * ERE * Z_k / lam(k);
		tmpmat.diag() += sig2_k;
		double logdet_val;
		double logdet_sign;
		log_det(logdet_val, logdet_sign, tmpmat);

		loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(tmpmat, resid_k));
	}
	loglik += 0.5 * (static_cast<double>(nT - 1 - std::abs(index2 - index1))) *
			  std::log(1.0 - std::pow((std::exp(2.0*zprho)-1.0)/(std::exp(2.0*zprho)+1.0),2.0))+
     	      2.0*zprho - 2.0*std::log(std::exp(2.0*zprho)+1.0);
	return loglik;
}
