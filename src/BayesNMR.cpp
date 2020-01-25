#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "Listbuilder.h"
#include "misc_nmr.h"
// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

// [[Rcpp::export]]
Rcpp::List BayesNMR(const arma::vec& y,
					const arma::vec& sd,
					const arma::mat& x,
					const arma::mat& z,
					const arma::uvec& ids,
					const arma::uvec& iarm,
					const arma::vec& npt,
					const double& nu, // = degrees of freedom
					const double& c01_inv,
					const double& c02_inv,
					const int& K,
					const int& nT,
					const int& ndiscard,
					const int& nskip,
					const int& nkeep) {
	using namespace arma;
	using namespace Rcpp;
	using namespace R;
	using namespace std;

	const int ns = y.n_elem;
	const int nx = x.n_cols;
	const int nz = z.n_cols;
	const vec sd2 = arma::pow(sd, 2.0);

	ivec narms(ns, fill::zeros);
	for (int i = 0; i < ns; ++i) {
		++narms(ids(i));
	}

	/*************************
	Parameters for adaptive MH
	*************************/
	const double obj_rate = 0.44;
	const int batch_length = std::min(50, ndiscard+nkeep*nskip);
	const int batch_total = (ndiscard + nskip*nkeep) / batch_length;

	vec eta_accepts(K, fill::zeros);
	mat eta_rates(K, batch_total, fill::zeros);
	vec eta_tunings(K);
	eta_tunings.fill(0.1);
	int batch_num = 0;

	vec phi_accepts(nz, fill::zeros);
	mat phi_rates(nz, batch_total, fill::zeros);
	vec phi_tunings(nz);
	phi_tunings.fill(0.1);

	mat Rho_accepts(nT,nT, fill::zeros);
	cube Rho_rates(nT,nT,batch_total, fill::zeros);
	mat Rho_tunings(nT,nT);
	Rho_tunings.fill(0.1);

	
	/********************
	Initialize parameters 
	********************/
	vec beta(nx+nT, fill::zeros);
	vec xb(ns, fill::zeros);
	vec resid = y - xb;
	vec Rgam(ns, fill::zeros);
	mat Rho(nT,nT, fill::zeros);
	mat pRho(nT,nT, fill::zeros);
	for (int i = 0; i < nT; ++i) {
		Rho(i,i) = 1.0;
		pRho(i,i) = 1.0;
	}
	vec sig2(ns, fill::ones);

	vec phi(nz, fill::ones);
	vec Z = arma::exp(z * phi);

	vec lam(K, fill::ones);


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




	/*******************
	Begin burn-in period
	*******************/
	Rcout << "Warming up" << endl;
	{
		Progress prog(ndiscard, true);
		for (int idiscard = 0; idiscard < ndiscard; ++idiscard) {
			R_CheckUserInterrupt();
			/**********
			Sample sig2
			**********/
			vec shape = 0.5 * npt + 0.00001;
			vec rate = 0.5 * npt % arma::pow(resid - Z % Rgam, 2.0) + 0.5 * (npt - 1.0) % sd2 + 0.00001;
			for (int i = 0; i < ns; ++i) {
				sig2(i) = rate(i) / ::Rf_rgamma(shape(i), 1.0);
			}

			/**********
			Sample beta
			**********/
			mat SigBetainv(nx+nT, nx+nT, fill::zeros);
			vec muBetaTmp(nx+nT, fill::zeros);
			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				mat Sig_k = arma::diagmat(sig2(idx) / npt(idx));
				vec y_k = y(idx);
				mat E_k = Eks(k);
				mat X_k = arma::join_horiz(Xks(k), E_k.t());
				mat Z_k = arma::diagmat(Z(idx));

				mat tmpmat = Z_k * E_k.t() * Rho * E_k * Z_k / lam(k) + Sig_k;
				SigBetainv += X_k.t() * arma::solve(tmpmat, X_k);
				muBetaTmp += X_k.t() * arma::solve(tmpmat, y_k);
			}
			SigBetainv.diag() += c01_inv;
			SigBetainv = 0.5 * (SigBetainv + SigBetainv.t());
			mat SigBetainvChol = chol(SigBetainv);
			vec muBeta = arma::solve(arma::trimatu(SigBetainvChol), arma::solve(trimatl(SigBetainvChol.t()), muBetaTmp));
			vec btmp(nx+nT);
			std::generate(btmp.begin(), btmp.end(), ::norm_rand);
			beta = muBeta + arma::solve(arma::trimatu(SigBetainvChol), btmp);
			for (int k=0; k < K; ++k) {
				mat E_k = Eks(k);
				mat X_k = arma::join_horiz(Xks(k), E_k.t());
				uvec idx = idxks(k);
				xb(idx) = X_k * beta;
			}
			resid = y - xb;

			/***********************
			Sample eta = log(lambda)
			***********************/
			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				vec sig2_k = sig2(idx) / npt(idx);
				vec resid_k = resid(idx);
				mat E_k = Eks(k);
				vec z_k = Z(idx);
				mat ERE_k = E_k.t() * Rho * E_k;

				double eta_k = std::log(lam(k));
				double eta_prop = ::norm_rand() * std::exp(eta_tunings(k)) + eta_k;
				
				// log-likelihood difference
				double ll_diff = loglik_eta(eta_prop, nu, resid_k, z_k, ERE_k, sig2_k) - 
								loglik_eta(eta_k, nu, resid_k, z_k, ERE_k, sig2_k);

				if (std::log(::unif_rand()) < ll_diff) {
					lam(k) = std::exp(eta_prop);
					++eta_accepts(k);
				}
			}


			/*********
			Sample phi
			*********/
			for (int g=0; g < nz; ++g) {
				vec phi_prop = phi;
				phi_prop(g) = ::norm_rand() * std::exp(phi_tunings(g)) + phi(g);

				// log-likelihood difference
				double ll_diff = loglik_phi(phi_prop, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks) -
						  loglik_phi(phi, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks);
				if (std::log(::unif_rand()) < ll_diff) {
					phi(g) = phi_prop(g);
					++phi_accepts(g);
				}
			}
			Z = arma::exp(z * phi);


			/*********
			Sample Rho
			*********/
			for (int iR=0; iR < nT-1; ++iR) {
				for (int iC=iR+1; iC < nT; ++iC) {
					double zprho = 0.5 * std::log((1.0 + pRho(iR,iC)) / (1.0 - pRho(iR,iC)));
					double zprho_prop = ::norm_rand() * std::exp(Rho_tunings(iR,iC)) + zprho;

					// log-likelihood difference
					double ll_diff = loglik_z(zprho_prop, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks) -
									 loglik_z(zprho, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks);
					if (std::log(::unif_rand()) < ll_diff) {
						pRho(iR,iC) = (std::exp(2.0 * zprho_prop) - 1.0) / (std::exp(2.0 * zprho_prop) + 1.0);
						pRho(iC,iR) = pRho(iR,iC);
						Rho = pRho_to_Rho(pRho);
						++Rho_accepts(iR,iC);
					}
				}
			}

			/***********
			Sample Rgam
			***********/
			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				mat Z_k = arma::diagmat(Z(idx));
				vec sig2_k = sig2(idx) / npt(idx);
				mat E_k = Eks(k);
				mat ERE = E_k.t() * Rho * E_k;
				vec resid_k = resid(idx);

				mat SigRgami = Z_k * diagmat(1.0 / sig2_k) * Z_k + lam(k) * ERE.i();
				SigRgami = 0.5 * (SigRgami + SigRgami.t());
				
				mat SigRgamiChol = chol(SigRgami);
				vec muRgam = arma::solve(arma::trimatu(SigRgamiChol), arma::solve(arma::trimatl(SigRgamiChol.t()), Z_k * diagmat(1.0 / sig2_k) * resid_k));

				int n_k = idx.n_elem;
				vec gtmp(n_k);
				std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
				Rgam(idx) = muRgam + arma::solve(arma::trimatu(SigRgamiChol), gtmp);
			}

			if ((idiscard+1) % batch_length == 0) {
				Rho_accepts /= static_cast<double>(batch_length);
				Rho_rates.slice(batch_num) = Rho_accepts;
				for (int iR=0; iR < nT-1; ++iR) {
					for (int iC=iR+1; iC < nT; ++iC) {
						if (Rho_accepts(iR,iC) > obj_rate) {
							Rho_tunings(iR,iC) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						} else {
							Rho_tunings(iR,iC) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						}
					}
				}

				eta_accepts /= static_cast<double>(batch_length);
				eta_rates.col(batch_num) = eta_accepts;
				for (int i=0; i<K; ++i) {
					if (eta_accepts(i) > obj_rate) {
						eta_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					} else {
						eta_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));

					}
				}

				phi_accepts /= static_cast<double>(batch_length);
				phi_rates.col(batch_num) = phi_accepts;
				for (int i=0; i<nz; ++i) {
					if (phi_accepts(i) > obj_rate) {
						phi_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					} else {
						phi_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
					}
				}
				++batch_num;
				Rho_accepts.fill(0.0);
				eta_accepts.fill(0.0);
				phi_accepts.fill(0.0);
			}
			prog.increment();
		}
	}

	/***********************
	Begin posterior sampling
	***********************/
	mat beta_save(nx+nT, nkeep);
	mat sig2_save(ns, nkeep);
	mat lam_save(K, nkeep);
	mat phi_save(nz, nkeep);
	cube Rho_save(nT,nT,nkeep);
	mat gam_save(ns, nkeep);


	int icount_mh = ndiscard;
	Rcout << "Saving posterior samples" << endl;
	{
		Progress prog(nkeep, true);
		for (int ikeep=0; ikeep < nkeep; ++ikeep) {
			R_CheckUserInterrupt();
			for (int iskip=0; iskip < nskip; ++iskip) {
				++icount_mh;
				/**********
				Sample sig2
				**********/
				vec shape = 0.5 * npt + 0.00001;
				vec rate = 0.5 * npt % arma::pow(resid - Z % Rgam, 2.0) + 0.5 * (npt - 1.0) % sd2 + 0.00001;
				for (int i = 0; i < ns; ++i) {
					sig2(i) = rate(i) / ::Rf_rgamma(shape(i), 1.0);
				}

				/**********
				Sample beta
				**********/
				mat SigBetainv(nx+nT, nx+nT, fill::zeros);
				vec muBetaTmp(nx+nT, fill::zeros);
				for (int k=0; k < K; ++k) {
					uvec idx = idxks(k);
					mat Sig_k = arma::diagmat(sig2(idx) / npt(idx));
					vec y_k = y(idx);
					mat E_k = Eks(k);
					mat X_k = arma::join_horiz(Xks(k), E_k.t());
					mat Z_k = arma::diagmat(Z(idx));

					mat tmpmat = Z_k * E_k.t() * Rho * E_k * Z_k / lam(k) + Sig_k;
					SigBetainv += X_k.t() * arma::solve(tmpmat, X_k);
					muBetaTmp += X_k.t() * arma::solve(tmpmat, y_k);
				}
				SigBetainv.diag() += c01_inv;
				SigBetainv = 0.5 * (SigBetainv + SigBetainv.t());
				mat SigBetainvChol = chol(SigBetainv);
				vec muBeta = arma::solve(arma::trimatu(SigBetainvChol), arma::solve(trimatl(SigBetainvChol.t()), muBetaTmp));
				vec btmp(nx+nT);
				std::generate(btmp.begin(), btmp.end(), ::norm_rand);
				beta = muBeta + arma::solve(arma::trimatu(SigBetainvChol), btmp);
				for (int k=0; k < K; ++k) {
					mat E_k = Eks(k);
					mat X_k = arma::join_horiz(Xks(k), E_k.t());
					uvec idx = idxks(k);
					xb(idx) = X_k * beta;
				}
				resid = y - xb;

				/***********************
				Sample eta = log(lambda)
				***********************/
				for (int k=0; k < K; ++k) {
					uvec idx = idxks(k);
					vec sig2_k = sig2(idx) / npt(idx);
					vec resid_k = resid(idx);
					mat E_k = Eks(k);
					vec z_k = Z(idx);
					mat ERE_k = E_k.t() * Rho * E_k;

					double eta_k = std::log(lam(k));
					double eta_prop = ::norm_rand() * std::exp(eta_tunings(k)) + eta_k;
					
					// log-likelihood difference
					double ll_diff = loglik_eta(eta_prop, nu, resid_k, z_k, ERE_k, sig2_k) - 
									loglik_eta(eta_k, nu, resid_k, z_k, ERE_k, sig2_k);

					if (std::log(::unif_rand()) < ll_diff) {
						lam(k) = std::exp(eta_prop);
						++eta_accepts(k);
					}
				}


				/*********
				Sample phi
				*********/
				for (int g=0; g < nz; ++g) {
					vec phi_prop = phi;
					phi_prop(g) = ::norm_rand() * std::exp(phi_tunings(g)) + phi(g);

					// log-likelihood difference
					double ll_diff = loglik_phi(phi_prop, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks) -
							  loglik_phi(phi, z, c02_inv, lam, sig2 / npt, Rho, resid, Eks, idxks);
					if (std::log(::unif_rand()) < ll_diff) {
						phi(g) = phi_prop(g);
						++phi_accepts(g);
					}
				}
				Z = arma::exp(z * phi);


				/*********
				Sample Rho
				*********/
				for (int iR=0; iR < nT-1; ++iR) {
					for (int iC=iR+1; iC < nT; ++iC) {
						double zprho = 0.5 * std::log((1.0 + pRho(iR,iC)) / (1.0 - pRho(iR,iC)));
						double zprho_prop = ::norm_rand() * std::exp(Rho_tunings(iR,iC)) + zprho;

						// log-likelihood difference
						double ll_diff = loglik_z(zprho_prop, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks) -
										 loglik_z(zprho, iR, iC, pRho, lam, sig2 / npt, Z, resid, Eks, idxks);
						if (std::log(::unif_rand()) < ll_diff) {
							pRho(iR,iC) = (std::exp(2.0 * zprho_prop) - 1.0) / (std::exp(2.0 * zprho_prop) + 1.0);
							pRho(iC,iR) = pRho(iR,iC);
							Rho = pRho_to_Rho(pRho);
							++Rho_accepts(iR,iC);
						}
					}
				}

				/***********
				Sample Rgam
				***********/
				for (int k=0; k < K; ++k) {
					uvec idx = idxks(k);
					mat Z_k = arma::diagmat(Z(idx));
					vec sig2_k = sig2(idx) / npt(idx);
					mat E_k = Eks(k);
					mat ERE = E_k.t() * Rho * E_k;
					vec resid_k = resid(idx);

					mat SigRgami = Z_k * diagmat(1.0 / sig2_k) * Z_k + lam(k) * ERE.i();
					SigRgami = 0.5 * (SigRgami + SigRgami.t());
					mat SigRgamiChol = cholmod(SigRgami);
					vec muRgam = arma::solve(arma::trimatu(SigRgamiChol), arma::solve(arma::trimatl(SigRgamiChol.t()), Z_k * diagmat(1.0 / sig2_k) * resid_k));

					int n_k = idx.n_elem;
					vec gtmp(n_k);
					std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
					Rgam(idx) = muRgam + arma::solve(arma::trimatu(SigRgamiChol), gtmp);
				}


				if (icount_mh % batch_length == 0) {
					Rho_accepts /= static_cast<double>(batch_length);
					Rho_rates.slice(batch_num) = Rho_accepts;
					for (int iR=0; iR < nT-1; ++iR) {
						for (int iC=iR+1; iC < nT; ++iC) {
							if (Rho_accepts(iR,iC) > obj_rate) {
								Rho_tunings(iR,iC) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
							} else {
								Rho_tunings(iR,iC) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
							}
						}
					}

					eta_accepts /= static_cast<double>(batch_length);
					eta_rates.col(batch_num) = eta_accepts;
					for (int i=0; i<K; ++i) {
						if (eta_accepts(i) > obj_rate) {
							eta_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						} else {
							eta_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));

						}
					}

					phi_accepts /= static_cast<double>(batch_length);
					phi_rates.col(batch_num) = phi_accepts;
					for (int i=0; i<nz; ++i) {
						if (phi_accepts(i) > obj_rate) {
							phi_tunings(i) += std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						} else {
							phi_tunings(i) -= std::min(0.01, 1.0 / std::sqrt(static_cast<double>(batch_num+1)));
						}
					}
					++batch_num;
					Rho_accepts.fill(0.0);
					eta_accepts.fill(0.0);
					phi_accepts.fill(0.0);
				}
			}
			beta_save.col(ikeep) = beta;
			phi_save.col(ikeep) = phi;
			lam_save.col(ikeep) = lam;
			sig2_save.col(ikeep) = sig2;
			Rho_save.slice(ikeep) = Rho;
			gam_save.col(ikeep) = Rgam;

			vec LL(K);
			for (int k=0; k < K; ++k) {
				uvec idx = idxks(k);
				int Tk = idx.n_elem;
				LL(k) = int_ObservedLik(Tk, resid(idx), lam(k), Z(idx), Eks(k), sig2(idx), npt(idx), Rho, nu);
			}
			double m = max(LL);
			

			prog.increment();
		}
	}

	return ListBuilder()
	.add("beta", beta_save)
	.add("phi", phi_save)
	.add("lam", lam_save)
	.add("sig2", sig2_save)
	.add("Rho", Rho_save)
	.add("gam", gam_save)
	.add("Rho_mh", Rho_rates)
	.add("eta_mh", eta_rates)
	.add("phi_mh", phi_rates);
}
