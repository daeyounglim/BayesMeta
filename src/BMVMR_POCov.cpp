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
#include "nelmin.h"
#include "ListBuilder.h"

// [[Rcpp::depends(RcppArmadillo,RcppProgress))]]

/*
// #include "misc_nmr.h"
C++ code for fitting
Bayesian Inference for Multivariate Meta-Regression With a Partially Observed Within-Study Sample Covariance Matrix

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
					   const arma::uvec& Treat,
					   const arma::uvec& Trial,
					   const arma::vec& Npt,
					   const double& c0,
					   const double& dj0, // hyperparameter for Omega
					   const double& d0, // hyperparameter for Sigma
					   const double& s0, // hyperparameter for Sigma
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
	mat pRtk(N, J * (J - 1) / 2, fill::zeros);
	for (int k = 0; k < K; ++k) {
		for (int t = 0; t < T; ++t) {
			Sigmainvs(t,k) = eye<mat>(J,J);
		}
	}
	mat delta(N, J, fill::ones);
	mat Rho(J, J, fill::eye);
	mat pRho(J, J, fill::eye);

	/************
	Miscellaneous
	************/
	const mat Omega0inv = arma::inv(Omega0);
	const mat Sigma0inv = arma::inv(Sigma0);
	const double df = s0 + arma::accu(Npt);
	const double shape_omega = static_cast<double>(K) + dj0;
	mat resid = Outcome;
	mat pR_rates(arma::size(pRtk), fill::zeros);
	mat pRho_rates(arma::size(pRho), fill::zeros);
	mat delta_rates(arma::size(delta), fill::zeros);

	/*********
	Containers
	*********/
	mat theta_save(nt, nkeep, fill::zeros);
	cube Omegainv_save(nw*J, nw*J, nkeep, fill::zeros);
	cube Sigmainv_save(J, J, nkeep, fill::zeros);
	cube delta_save(N, J, nkeep, fill::zeros);
	cube Rho_save(J, J, nkeep, fill::zeros);
	cube Rtk_save(N, J * (J - 1) / 2, nkeep, fill::zeros);

	/*******************
	Begin burn-in period
	*******************/
	int icount = 0;
	if (verbose) {
		Rcout << "Warming up" << endl;
	}
	{
		Progress prog(ndiscard, verbose);
		for (int idiscard = 0; idiscard < ndiscard+nkeep; ++idiscard) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}


			/***********
			Update theta
			***********/
			mat Sig_theta(nt, nt, fill::zeros);
			Sig_theta.diag().fill(1.0 / c0);
			vec mu_theta(nt, fill::zeros);
			for (int k = 0; k < K; ++k) {
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				mat XSX(nt, nt, fill::zeros);
				mat WSX(nw*J, nt, fill::zeros);
				vec WSy(nw*J, fill::zeros);
				vec XSy(nt, fill::zeros);
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
					XSX += ntk * X.t() * Sigmainv * X;
					XSy += ntk * X.t() * Sigmainv * y_i.t();
					mat WS = ntk * W.t() * Sigmainv;
					Sig_gamk += WS * W;
					WSX += WS * X;
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
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = XCovariate.row(i_k);
					rowvec w_i = WCovariate.row(i_k);
					rowvec y_i = Outcome.row(i_k);
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
				mat qq = Omega0inv;
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
						double ntk = Npt(i);
						vec gam_k = gamR.col(k);
						mat R = veclinv(trans(pRtk.row(i)), J);
						R.diag().fill(1.0);
						R = pRho_to_Rho(R);
						mat V = arma::diagmat(SD.row(i));
						rowvec w_i = WCovariate.row(i);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
						qq += ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;
					}
					mat Siginv_new = rwish(df, qq);
					for (int t = 0; t < T; ++t) {
						for (int k = 0; k < K; ++k) {
							Sigmainvs(t, k) = Siginv_new;
						}
					}
				} else if (fmodel == 3) {
					// Update sig2
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						double ntk = Npt(i);
						rowvec w_i = WCovariate.row(i);
						mat R = veclinv(trans(pRtk.row(i)), J);
						R.diag().fill(1.0);
						R = pRho_to_Rho(R);
						vec gam_k = gamR.col(k);
						mat V = arma::diagmat(SD.row(i));
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
						mat qq = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;

						rowvec delta_i = delta.row(i);

						for (int j = 0; j < J; ++j) {
							auto fx_delta = [&](double delta_input[])->double {
								double deltap = delta_input[0];
								double sig2p = std::exp(deltap);
								rowvec sig2_ip = delta_i;
								sig2_ip(j) = sig2p;
								mat Vinvp = arma::diagmat(1.0 / sig2_ip);
								mat Sigmainvp = Vinvp * arma::inv_sympd(Rho) * Vinvp;

								double loglik = -0.5 * dot(qq, Sigmainvp) - (s0 + ntk) * deltap - d0 * std::exp(-deltap);
							    return -loglik;
							};

							double dstar = std::log(delta_i(j));
							double start[] = { dstar };
							double xmin[] = { 0.0 };
							double ynewlo = 0.0;
							double reqmin = 1.0e-20;
							int konvge = 5;
							int kcount = 1000;
							double step[] = { 0.3 };
							int icount = 0;
							int numres = 0;
							int ifault = 0;
							nelmin(fx_delta, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
							double xmax = xmin[0];
							double minll = ynewlo;

							mat cl(5,3, fill::zeros);
							vec dl(5, fill::zeros);
							double step_size = 0.2;
							delta_burnin_block:
								for (int iii=0; iii < 5; ++iii) {
									double e1 = static_cast<double>(iii-2);
									cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
									cl(iii,1) = xmax + e1 * step_size;
									cl(iii,2) = 1.0;
									start[0] = xmax + e1 * step_size;
									dl(iii) = fx_delta(start);
								}

							for (int ni=0; ni < 5; ++ni) {
								if ((ni+1) != 3) {
									if (dl(ni) <= minll) {
										step_size *= 1.2;
										goto delta_burnin_block;
									}
								}
							}

							vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
							double sigmaa = std::sqrt(0.5 / fl(0));


							double dstar_prop = ::norm_rand() * sigmaa + xmax;
							// log-likelihood difference
							start[0] = dstar;
							double ll_diff = fx_delta(start);
							start[0] = dstar_prop;
							ll_diff += -fx_delta(start)
									    - 0.5 * (std::pow(dstar - xmax, 2.0) - std::pow(dstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
							if (std::log(::unif_rand()) < ll_diff) {
								delta_i(j) = std::exp(dstar_prop);
								delta(i, j) = delta_i(j);
								++delta_rates(i,j);
							}
						}
					}

					// Update Rho
					mat qq(J, J, fill::zeros);
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						double ntk = Npt(i);
						rowvec w_i = WCovariate.row(i);
						mat R = veclinv(trans(pRtk.row(i)), J);
						R.diag().fill(1.0);
						R = pRho_to_Rho(R);
						vec gam_k = gamR.col(k);
						mat V = arma::diagmat(SD.row(i));
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							W(j, span(j*nw, (j+1)*nw-1)) = w_i;
						}
						vec resid_i = arma::trans(resid.row(i)) - W * gam_k;
						mat dAd = ntk * resid_i * resid_i.t() + (ntk - 1.0) * V * R * V;
						mat siginvm = arma::diagmat(1.0 / delta.row(i));
						qq += siginvm * dAd * siginvm;
					}
					for (int iR = 0; iR < J-1; ++iR) {
						for (int iC =iR+1; iC < J; ++iC) {
							double zprho = 0.5 * std::log((1.0 + pRho(iR,iC)) / (1.0 - pRho(iR,iC)));
							auto fx_zrho = [&](double zprho_input[])->double {
								double zprhop = zprho_input[0];
								mat pRhop = pRho;
								pRhop(iR, iC) = (std::exp(2.0 * zprhop) - 1.0) / (std::exp(2.0 * zprhop) + 1.0);
								pRhop(iC, iR) = pRhop(iR, iC);
								mat Rhop = pRho_to_Rho(pRhop);
								double logdet_val;
								double logdet_sign;
								log_det(logdet_val, logdet_sign, Rhop);

								double loglik = -0.5 * dot(qq, arma::inv_sympd(Rhop)) - 0.5 * df * logdet_val;
								loglik += 0.5 * (static_cast<double>(J - 1 - std::abs(iC - iR))) *
										  std::log(1.0 - std::pow((std::exp(2.0*zprhop)-1.0)/(std::exp(2.0*zprhop)+1.0),2.0))+
							     	      2.0*zprhop - 2.0*std::log(std::exp(2.0*zprhop)+1.0);
							    return -loglik;
							};
							double start[] = { zprho };
							double xmin[] = { 0.0 };
							double ynewlo = 0.0;
							double reqmin = 1.0e-20;
							int konvge = 5;
							int kcount = 1000;
							double step[] = { 0.03 };
							int icount = 0;
							int numres = 0;
							int ifault = 0;
							nelmin(fx_zrho, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
							double xmax = xmin[0];
							double minll = ynewlo;

							mat cl(5,3, fill::zeros);
							vec dl(5, fill::zeros);
							double step_size = 0.2;
							pRho_burnin_block:
								for (int iii=0; iii < 5; ++iii) {
									double e1 = static_cast<double>(iii-2);
									cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
									cl(iii,1) = xmax + e1 * step_size;
									cl(iii,2) = 1.0;
									start[0] = xmax + e1 * step_size;
									dl(iii) = fx_zrho(start);
								}

							for (int ni=0; ni < 5; ++ni) {
								if ((ni+1) != 3) {
									if (dl(ni) <= minll) {
										step_size *= 1.2;
										goto pRho_burnin_block;
									}
								}
							}

							vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
							double sigmaa = std::sqrt(0.5 / fl(0));


							double zprho_prop = ::norm_rand() * sigmaa + xmax;
							// log-likelihood difference
							start[0] = zprho;
							double ll_diff = fx_zrho(start);
							start[0] = zprho_prop;
							ll_diff += -fx_zrho(start)
									    - 0.5 * (std::pow(zprho - xmax, 2.0) - std::pow(zprho_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
							if (std::log(::unif_rand()) < ll_diff) {
								pRho(iR,iC) = (std::exp(2.0 * zprho_prop) - 1.0) / (std::exp(2.0 * zprho_prop) + 1.0);
								pRho(iC,iR) = pRho(iR,iC);
								Rho = pRho_to_Rho(pRho);
								++pRho_rates(iR,iC);
							}
						}
					}

					// Update Sigmainvs
					mat Rhoinv = arma::inv_sympd(Rho);
					for (int i = 0; i < N; ++i) {
						int k = Trial(i);
						int t = Treat(i);
						mat siginvm = arma::diagmat(1.0 / delta.row(i));
						Sigmainvs(t, k) = siginvm * Rhoinv * siginvm;
					}

				} 
				// else if (fmodel == 4) {
				//
				// }

				// Update R
				for (int i = 0; i < N; ++i) {
					int k = Trial(i);
					int t = Treat(i);

					rowvec y_i = Outcome.row(i);
					rowvec prtk = pRtk.row(i);
					mat R = veclinv(trans(pRtk.row(i)), J);
					R.diag().fill(1.0);
					R = pRho_to_Rho(R);
					mat V = arma::diagmat(SD.row(i));
					double ntk = Npt(i);
					mat Sigmainv = Sigmainvs(t, k);

					for (int k = 0; k < J*(J-1)/2; ++k) {
						auto fx_rstar = [&](double rstar_input[])->double {
							double rstar = rstar_input[0];
							double zprho = (std::exp(2.0 * rstar) - 1.0) / (std::exp(2.0 * rstar) + 1.0);
							rowvec tmp_prtk = prtk;
							tmp_prtk(k) = zprho;
							mat R_prop = veclinv(trans(tmp_prtk), J);
							R_prop.diag().fill(1.0);
							R_prop = pRho_to_Rho(R_prop);
							double logdet_val;
							double logdet_sign;
							log_det(logdet_val, logdet_sign, R_prop);
							int iR = J - 2 - static_cast<int>(std::sqrt(-8.0*static_cast<double>(k) + 4.0*static_cast<double>(J*(J-1)-7))/2.0 - 0.5); // row index
							int iC = k + iR + 1 - J*(J-1)/2 + (J-iR)*((J-iR)-1)/2; // column index

							mat S_prop = V * R_prop * V;
							double loglik = -0.5 * (ntk - static_cast<double>(J) - 2.0) * logdet_val - 0.5 * (ntk - 1.0) * dot(S_prop, Sigmainv);
							loglik += 0.5 * (static_cast<double>(J - 1 - std::abs(iC - iR))) *
									  std::log(1.0 - std::pow((std::exp(2.0*zprho)-1.0)/(std::exp(2.0*zprho)+1.0),2.0))+
						     	      2.0*zprho - 2.0*std::log(std::exp(2.0*zprho)+1.0);
						    return -loglik;
						};
						double zstar = 0.5 * std::log((1.0 + prtk(k)) / (1.0 - prtk(k)));
						double start[] = { zstar };
						double xmin[] = { 0.0 };
						double ynewlo = 0.0;
						double reqmin = 1.0e-20;
						int konvge = 5;
						int kcount = 1000;
						double step[] = { 0.03 };
						int icount = 0;
						int numres = 0;
						int ifault = 0;
						nelmin(fx_rstar, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);					
						double xmax = xmin[0];
						double minll = ynewlo;

						mat cl(5,3, fill::zeros);
						vec dl(5, fill::zeros);
						double step_size = 0.2;
						R_burnin_block:
							for (int iii=0; iii < 5; ++iii) {
								double e1 = static_cast<double>(iii-2);
								cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
								cl(iii,1) = xmax + e1 * step_size;
								cl(iii,2) = 1.0;
								start[0] = xmax + e1 * step_size;
								dl(iii) = fx_rstar(start);
							}

						for (int ni=0; ni < 5; ++ni) {
							if ((ni+1) != 3) {
								if (dl(ni) <= minll) {
									step_size *= 1.2;
									goto R_burnin_block;
								}
							}
						}

						vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
						double sigmaa = std::sqrt(0.5 / fl(0));


						double zstar_prop = ::norm_rand() * sigmaa + xmax;
						// log-likelihood difference
						start[0] = zstar;
						double ll_diff = fx_rstar(start);
						start[0] = zstar_prop;
						ll_diff += -fx_rstar(start)
								    - 0.5 * (std::pow(zstar - xmax, 2.0) - std::pow(zstar_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
						if (std::log(::unif_rand()) < ll_diff) {
							prtk(k) = (std::exp(2.0 * zstar_prop) - 1.0) / (std::exp(2.0 * zstar_prop) + 1.0);
							pRtk.row(i) = prtk;
							++pR_rates(i,k);
						}
					}
				}
			}
			if (idiscard > ndiscard) {
				theta_save.col(icount) = theta;
				Omegainv_save.slice(icount) = Omegainv;

				if (fmodel == 1) {
					Sigmainv_save.slice(icount) = Sigmainvs(0,0);
				} else {
					if (fmodel == 2) {
						Sigmainv_save.slice(icount) = Sigmainvs(0,0);
					} else if (fmodel == 3) {
						delta_save.slice(icount) = delta;
						Rho_save.slice(icount) = Rho;
					}

					mat Rtk(arma::size(pRtk), fill::zeros);
					for (int i = 0; i < N; ++i) {
						mat R = veclinv(trans(pRtk.row(i)), J);
						R.diag().fill(1.0);
						R = pRho_to_Rho(R);
						Rtk.row(i) = arma::trans(vecl(R));
					}
					Rtk_save.slice(icount) = Rtk;
				}
				++icount;
			}
		}
	}

	ListBuilder out;
	if (fmodel == 1) {
		out = ListBuilder()
			.add("theta", theta_save)
			.add("Omegainv", Omegainv_save)
			.add("Sigmainv", Sigmainv_save);
	} else {
		if (fmodel == 2) {
			out = ListBuilder()
				.add("theta", theta_save)
				.add("Omegainv", Omegainv_save)
			    .add("Sigmainv", Sigmainv_save)
			    .add("R", Rtk_save)
			    .add("pR_acceptance", pR_rates);
		} else if (fmodel == 3) {
			out = ListBuilder()
				.add("theta", theta_save)
				.add("Omegainv", Omegainv_save)
				.add("delta", delta_save)
				.add("Rho", Rho_save)
				.add("delta_acceptance", delta_rates)
				.add("pRho_acceptance", pRho_rates)
				.add("R", Rtk_save)
			    .add("pR_acceptance", pR_rates);
		} else {
			out = ListBuilder().add("warning", "Returned without a result. Please check your fmodel.");
		}
	}
	return out;
}