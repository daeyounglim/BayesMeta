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
#include "rtmvnCpp.h"
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
Rcpp::List BMVSMR(const arma::mat& Outcome,
					const arma::mat& Covariate,
					const arma::vec& Treat,
					const arma::uvec& Patient,
					const arma::uvec& Trial,
					const double& a0,
					const double& b0,
					const double& c1,
					const double& c2,
					const double& c3,
					const double& d0,
					const arma::mat& S0,
					const double& d1,
					const arma::mat& S1,
					const double& v0,
					const double& tau,
					const int& K,
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
	const int xcols = Covariate.n_cols;
	const int nw = 2;
	const int nt = (xcols + nw) * J;
	const int ntd = nt + J;

	arma::field<arma::uvec> idxks(K);
	for (int k = 0; k < K; ++k) {
		uvec idx = find(Trial == k);
		idxks(k) = idx;
	}

	vec theta(nt, fill::zeros);
	vec lambda(N, fill::ones);
	vec delta(J, fill::ones);
	vec phi(N, fill::ones); // `zeta` in FORTRAN code
	mat zz(N, J, fill::ones);
	vec Ez(J, fill::ones);
	mat Sigmainv(J, J, fill::eye);
	mat Omegainv(nw*J, nw*J, fill::eye);

	mat gamR(nw*J, K, fill::zeros); // `xi` in FORTRAN code
	double nu = v0 + 1.0; // `vv` in FORTRAN code
	vec Q0inv(ntd, fill::zeros);
	Q0inv.head(xcols*J).fill(1.0 / c1);
	for (int i = xcols*J; i < nt; ++i) {
		Q0inv(i) = 1.0 / c2;
	}
	Q0inv.tail(J).fill(1.0 / c3);

	double shape_omega = static_cast<double>(K) + d1;
	double shape_sigma = static_cast<double>(N) + d0;
	mat S0inv = arma::inv(S0);
	mat S1inv = arma::inv(S1);
	double nu_rates = 0.0;
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
			for (int k = 0; k < K; ++k) {
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				vec gamR_k = gamR.col(k);
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = Covariate.row(i_k);
					rowvec w_i(2, fill::ones);
					rowvec y_i = Outcome.row(i_k);
					w_i(1) = Treat(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*nw, (j+1)*nw-1)) = w_i;
					}
					mat Az = lambda(i_k) * DSD;
					mat Bz = lambda(i_k) * DS * (y_i.t() - arma::join_horiz(X, W) * theta - W * gamR_k + delta % Ez) - phi(i_k);
					Az = 0.5 * (Az + Az.t());
					// mat Azinv = Az.i();
					// vec muz = Azinv * Bz;
					for (int j1 = 0; j1 < J; ++j1) {
						double asigma = Az(j1,j1);
						double amean = Bz(j1);
						for (int j2 = 0; j2 < J; ++j2) {
							if (j2 != j1) {
								amean -= Az(j1, j2) * zz(i_k, j2);
							}
						}
						double zmean = amean / asigma;
						double zsigma = 1.0 / asigma;

						double trim = -zmean / std::sqrt(zsigma);
						double rv = ltnormrnd(0.0, 1.0, trim);
						zz(i_k, j1) = zmean + rv * std::sqrt(zsigma);
					}
				}
			}

			/*********
			Update phi
			*********/
			for (int i = 0; i < N; ++i) {
				rowvec z_i = zz.row(i);
				phi(i) = ::Rf_rgamma(tau + static_cast<double>(J+1), 1.0) / (tau + arma::accu(z_i));
			}

			/***********
			Update alpha
			***********/
			mat Sig_alpha(ntd, ntd, fill::zeros);
			Sig_alpha.diag() += Q0inv;
			vec mu_alpha(ntd, fill::zeros);
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
					rowvec x_i = Covariate.row(i_k);
					rowvec w_i(2, fill::ones);
					rowvec y_i = Outcome.row(i_k);
					w_i(1) = Treat(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					mat Xstar(J, ntd, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*2, (j+1)*2-1)) = w_i;
						Xstar(j, ntd - J + j) = zz(i_k,j)-1.0; // Ez = 1.0
					}
					Xstar.head_cols(xcols*J) = X;
					Xstar.cols(xcols*J, ntd - J - 1) = W;

					XSX += lambda(i_k) * Xstar.t() * Sigmainv * Xstar;
					XSy += lambda(i_k) * Xstar.t() * Sigmainv * y_i.t();
					mat WS = lambda(i_k) * W.t() * Sigmainv;
					Sig_gamk += WS * W;
					WSX += WS * Xstar;
					WSy += WS * y_i.t();
				}
				mat Sig_gamk_inv = arma::inv(Sig_gamk);
				Sig_alpha += XSX - WSX.t() * Sig_gamk_inv * WSX;
				mu_alpha += XSy - WSX.t() * Sig_gamk_inv * WSy;
			}
			Sig_alpha = 0.5 * (Sig_alpha + Sig_alpha.t());
			mat Sig_alphaChol = chol(Sig_alpha);
			vec alpha = arma::solve(arma::trimatu(Sig_alphaChol), arma::solve(trimatl(Sig_alphaChol.t()), mu_alpha));
			vec atmp(ntd);
			std::generate(atmp.begin(), atmp.end(), ::norm_rand);
			alpha += arma::solve(arma::trimatu(Sig_alphaChol), atmp);
			theta = alpha.head(nt);
			delta = alpha.tail(J);

			/**************
			Update gamR(xi)
			**************/
			for (int k = 0; k < K; ++k) {
				mat Siggam = Omegainv;
				vec mugam(nw*J, fill::zeros);
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = Covariate.row(i_k);
					rowvec w_i(2, fill::ones);
					rowvec y_i = Outcome.row(i_k);
					rowvec z_i = zz.row(i_k) - 1.0;
					w_i(1) = Treat(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*2, (j+1)*2-1)) = w_i;
					}
					mat WS = W.t() * Sigmainv;
					Siggam += lambda(i_k) * WS * W;
					mugam += lambda(i_k) * WS * (y_i.t() - arma::join_horiz(X,W) * theta - delta % z_i.t());
				}
				Siggam = 0.5 * (Siggam + Siggam.t());
				mat SiggamChol = chol(Siggam);
				vec gam_k = arma::solve(arma::trimatu(SiggamChol), arma::solve(trimatl(SiggamChol.t()), mugam));
				vec gtmp(arma::size(gam_k));
				std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
				gamR.col(k) = gam_k + arma::solve(arma::trimatu(SiggamChol), gtmp);
			}

			/***********
			Update Omega
			***********/
			for (int jj = 0; jj < J; ++jj) {
				mat qq = S1inv;
				vec gamstar(2);
				for (int k = 0; k < K; ++k) {
					for (int j = 0; j < 2; ++j) {
						gamstar(j) = gamR(nw*jj+j, k);
					}

					for (int j1 = 0; j1 < 2; ++j1) {
						for (int j2 = 0; j2 < 2; ++j2) {
							qq(j1, j2) += gamstar(j1) * gamstar(j2);
						}
					}
				}
				Omegainv(span(jj*2, (jj+1)*2-1), span(jj*2, (jj+1)*2-1)) = rwish(shape_omega, qq.i());
			}

			/***********
			Update Sigma
			***********/
			mat Asig(J, J, fill::zeros);
			for (int k = 0; k < K; ++k) {
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				vec gamR_k = gamR.col(k);
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = Covariate.row(i_k);
					rowvec w_i(2, fill::ones);
					rowvec y_i = Outcome.row(i_k);
					rowvec z_i = zz.row(i_k) - 1.0;
					w_i(1) = Treat(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*2, (j+1)*2-1)) = w_i;
					}
					vec resid = y_i.t() - arma::join_horiz(X,W) * theta - W * gamR_k - delta % z_i.t();
					Asig += lambda(i_k) * resid * resid.t();
				}
			}
			Sigmainv = rwish(shape_sigma, arma::inv(S0inv + Asig));

			/********
			Update nu
			********/
			auto fx_vv = [&](double star[])->double {
				double vv = v0 + std::exp(star[0]);
				double loglik = static_cast<double>(N) * (R::lgammafn(0.5 * (vv + static_cast<double>(J))) - R::lgammafn(0.5 * vv))
							 + (0.5 * static_cast<double>(N) * vv  + a0 - 1.0) * std::log(vv) - b0 * vv + star[0];
				for (int k = 0; k < K; ++k) {
					uvec idxk = idxks(k);
					int n_k = idxk.n_elem;
					vec gamR_k = gamR.col(k);
					for (int i = 0; i < n_k; ++i) {
						int i_k = idxk(i);
						rowvec x_i = Covariate.row(i_k);
						rowvec w_i(2, fill::ones);
						rowvec y_i = Outcome.row(i_k);
						rowvec z_i = zz.row(i_k) - 1.0;
						w_i(1) = Treat(i_k);
						mat X(J, xcols * J, fill::zeros);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
							W(j, span(j*2, (j+1)*2-1)) = w_i;
						}
						vec resid = y_i.t() - arma::join_horiz(X,W) * theta - W * gamR_k - delta % z_i.t();
						loglik -= 0.5 * (vv + static_cast<double>(J)) * std::log(vv + arma::dot(resid, Sigmainv * resid));
					}
				}
				return -loglik;
			};
			double bold = std::log(nu - v0);
			double start[] = { bold };
			double xmin[] = { 0.0 };
			double ynewlo = 0.0;
			double reqmin = 1.0e-10;
			int konvge = 5;
			int kcount = 1000;
			double step[] = { 0.2 };
			int icount = 0;
			int numres = 0;
			int ifault = 0;
			nelmin(fx_vv, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
			double xmax = xmin[0];
			double minll = ynewlo;

			mat cl(5,3, fill::zeros);
			vec dl(5, fill::zeros);
			double step_size = 0.1;
			vv_burnin_block:
				for (int iii=0; iii < 5; ++iii) {
					double e1 = static_cast<double>(iii-2);
					cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
					cl(iii,1) = xmax + e1 * step_size;
					cl(iii,2) = 1.0;
					start[0] = xmax + e1 * step_size;
					dl(iii) = fx_vv(start);
				}

			for (int ni=0; ni < 5; ++ni) {
				if ((ni+1) != 3) {
					if (dl(ni) <= minll) {
						step_size *= 1.2;
						goto vv_burnin_block;
					}
				}
			}

			vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
			double sigmaa = std::sqrt(0.5 / fl(0));


			double vv_prop = ::norm_rand() * sigmaa + xmax;
			start[0] = vv_prop;
			double ll_diff = fx_vv(start) - 0.5 * (std::pow(bold - xmax, 2.0) - std::pow(vv_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
			start[0] = bold;
			ll_diff -= fx_vv(start);
			if (std::log(::unif_rand()) < ll_diff) {
				nu = v0 + std::exp(vv_prop);
				++nu_rates;
			}

			/************
			Update lambda
			************/
			double shape_lam = 0.5 * (nu + static_cast<double>(J));
			for (int k = 0; k < K; ++k) {
				uvec idxk = idxks(k);
				int n_k = idxk.n_elem;
				vec gamR_k = gamR.col(k);
				for (int i = 0; i < n_k; ++i) {
					int i_k = idxk(i);
					rowvec x_i = Covariate.row(i_k);
					rowvec w_i(2, fill::ones);
					rowvec y_i = Outcome.row(i_k);
					w_i(1) = Treat(i_k);
					mat X(J, xcols * J, fill::zeros);
					mat W(J, nw*J, fill::zeros);
					rowvec z_i = zz.row(i_k) - 1.0;

					for (int j = 0; j < J; ++j) {
						X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
						W(j, span(j*2, (j+1)*2-1)) = w_i;
					}
					vec resid = y_i.t() - arma::join_horiz(X,W) * theta - W * gamR_k - delta % z_i.t();
					lambda(i_k) = ::Rf_rgamma(shape_lam, 1.0) / (0.5 * (nu + arma::dot(resid, resid)));
				}
			}
		}
	}

	mat thetaSave(nt, nkeep, fill::zeros);
	mat lambdaSave(N, nkeep, fill::zeros);
	mat deltaSave(J, nkeep, fill::zeros);
	vec nuSave(nkeep, fill::zeros);

	if (verbose) {
		Rcout << "Saving posterior samples" << endl;
	}
	{
		Progress prog(nkeep, verbose);
		for (int ikeep = 0; ikeep < nkeep; ++ikeep) {
			if (Progress::check_abort()) {
				return Rcpp::List::create(Rcpp::Named("error") = "user interrupt aborted");
			}
			for (int iskip = 0; iskip < nskip; ++iskip) {

				/********
				Update zz
				********/
				mat DS = arma::diagmat(delta) * Sigmainv;
				mat DSD = DS * arma::diagmat(delta);
				for (int k = 0; k < K; ++k) {
					uvec idxk = idxks(k);
					int n_k = idxk.n_elem;
					vec gamR_k = gamR.col(k);
					for (int i = 0; i < n_k; ++i) {
						int i_k = idxk(i);
						rowvec x_i = Covariate.row(i_k);
						rowvec w_i(2, fill::ones);
						rowvec y_i = Outcome.row(i_k);
						w_i(1) = Treat(i_k);
						mat X(J, xcols * J, fill::zeros);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
							W(j, span(j*2, (j+1)*2-1)) = w_i;
						}
						mat Az = lambda(i_k) * DSD;
						mat Bz = lambda(i_k) * DS * (y_i.t() - arma::join_horiz(X, W) * theta - W * gamR_k + delta % Ez) - phi(i_k);
						Az = 0.5 * (Az + Az.t());
						
						for (int j1 = 0; j1 < J; ++j1) {
							double asigma = Az(j1,j1);
							double amean = Bz(j1);
							for (int j2 = 0; j2 < J; ++j2) {
								if (j2 != j1) {
									amean -= Az(j1, j2) * zz(i_k, j2);
								}
							}
							double zmean = amean / asigma;
							double zsigma = 1.0 / asigma;

							double trim = -zmean / std::sqrt(zsigma);
							double rv = ltnormrnd(0.0, 1.0, trim);
							zz(i_k, j1) = zmean + rv * std::sqrt(zsigma);
						}
					}
				}

				/*********
				Update phi
				*********/
				for (int i = 0; i < N; ++i) {
					rowvec z_i = zz.row(i);
					phi(i) = ::Rf_rgamma(tau + static_cast<double>(J+1), 1.0) / (tau + arma::accu(z_i));
				}

				/***********
				Update alpha
				***********/
				mat Sig_alpha(ntd, ntd, fill::zeros);
				Sig_alpha.diag() += Q0inv;
				vec mu_alpha(ntd, fill::zeros);
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
						rowvec x_i = Covariate.row(i_k);
						rowvec w_i(2, fill::ones);
						rowvec y_i = Outcome.row(i_k);
						w_i(1) = Treat(i_k);
						mat X(J, xcols * J, fill::zeros);
						mat W(J, nw*J, fill::zeros);
						mat Xstar(J, ntd, fill::zeros);
						for (int j = 0; j < J; ++j) {
							X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
							W(j, span(j*2, (j+1)*2-1)) = w_i;
							Xstar(j, ntd - J + j) = zz(i_k,j)-1.0; // Ez = 1.0
						}
						Xstar.head_cols(xcols*J) = X;
						Xstar.cols(xcols*J, ntd - J - 1) = W;

						XSX += lambda(i_k) * Xstar.t() * Sigmainv * Xstar;
						XSy += lambda(i_k) * Xstar.t() * Sigmainv * y_i.t();
						mat WS = lambda(i_k) * W.t() * Sigmainv;
						Sig_gamk += WS * W;
						WSX += WS * Xstar;
						WSy += WS * y_i.t();
					}
					mat Sig_gamk_inv = arma::inv(Sig_gamk);
					Sig_alpha += XSX - WSX.t() * Sig_gamk_inv * WSX;
					mu_alpha += XSy - WSX.t() * Sig_gamk_inv * WSy;
				}
				Sig_alpha = 0.5 * (Sig_alpha + Sig_alpha.t());
				mat Sig_alphaChol = chol(Sig_alpha);
				vec alpha = arma::solve(arma::trimatu(Sig_alphaChol), arma::solve(trimatl(Sig_alphaChol.t()), mu_alpha));
				vec atmp(ntd);
				std::generate(atmp.begin(), atmp.end(), ::norm_rand);
				alpha += arma::solve(arma::trimatu(Sig_alphaChol), atmp);
				theta = alpha.head(nt);
				delta = alpha.tail(J);


				/**************
				Update gamR(xi)
				**************/
				for (int k = 0; k < K; ++k) {
					mat Siggam = Omegainv;
					vec mugam(nw*J, fill::zeros);
					uvec idxk = idxks(k);
					int n_k = idxk.n_elem;
					for (int i = 0; i < n_k; ++i) {
						int i_k = idxk(i);
						rowvec x_i = Covariate.row(i_k);
						rowvec w_i(2, fill::ones);
						rowvec y_i = Outcome.row(i_k);
						rowvec z_i = zz.row(i_k) - 1.0;
						w_i(1) = Treat(i_k);
						mat X(J, xcols * J, fill::zeros);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
							W(j, span(j*2, (j+1)*2-1)) = w_i;
						}
						mat WS = W.t() * Sigmainv;
						Siggam += lambda(i_k) * WS * W;
						mugam += lambda(i_k) * WS * (y_i.t() - arma::join_horiz(X,W) * theta - delta % z_i.t());
					}
					Siggam = 0.5 * (Siggam + Siggam.t());
					mat SiggamChol = chol(Siggam);
					vec gam_k = arma::solve(arma::trimatu(SiggamChol), arma::solve(trimatl(SiggamChol.t()), mugam));
					vec gtmp(arma::size(gam_k));
					std::generate(gtmp.begin(), gtmp.end(), ::norm_rand);
					gamR.col(k) = gam_k + arma::solve(arma::trimatu(SiggamChol), gtmp);
				}
				/***********
				Update Omega
				***********/
				for (int jj = 0; jj < J; ++jj) {
					mat qq = S1inv;
					vec gamstar(2);
					for (int k = 0; k < K; ++k) {
						for (int j = 0; j < 2; ++j) {
							gamstar(j) = gamR(nw*jj+j, k);
						}

						for (int j1 = 0; j1 < 2; ++j1) {
							for (int j2 = 0; j2 < 2; ++j2) {
								qq(j1, j2) += gamstar(j1) * gamstar(j2);
							}
						}
					}
					Omegainv(span(jj*2, (jj+1)*2-1), span(jj*2, (jj+1)*2-1)) = rwish(shape_omega, qq.i());
				}

				/***********
				Update Sigma
				***********/
				mat Asig(J, J, fill::zeros);
				for (int k = 0; k < K; ++k) {
					uvec idxk = idxks(k);
					int n_k = idxk.n_elem;
					vec gamR_k = gamR.col(k);
					for (int i = 0; i < n_k; ++i) {
						int i_k = idxk(i);
						rowvec x_i = Covariate.row(i_k);
						rowvec w_i(2, fill::ones);
						rowvec y_i = Outcome.row(i_k);
						rowvec z_i = zz.row(i_k) - 1.0;
						w_i(1) = Treat(i_k);
						mat X(J, xcols * J, fill::zeros);
						mat W(J, nw*J, fill::zeros);
						for (int j = 0; j < J; ++j) {
							X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
							W(j, span(j*2, (j+1)*2-1)) = w_i;
						}
						vec resid = y_i.t() - arma::join_horiz(X,W) * theta - W * gamR_k - delta % z_i.t();
						Asig += lambda(i_k) * resid * resid.t();
					}
				}
				Sigmainv = rwish(shape_sigma, arma::inv(S0inv + Asig));

				/********
				Update nu
				********/
				auto fx_vv = [&](double star[])->double {
					double vv = v0 + std::exp(star[0]);
					double loglik = static_cast<double>(N) * (R::lgammafn(0.5 * (vv + static_cast<double>(J))) - R::lgammafn(0.5 * vv))
								 + (0.5 * static_cast<double>(N) * vv  + a0 - 1.0) * std::log(vv) - b0 * vv + star[0];
					for (int k = 0; k < K; ++k) {
						uvec idxk = idxks(k);
						vec gamR_k = gamR.col(k);
						int n_k = idxk.n_elem;
						for (int i = 0; i < n_k; ++i) {
							int i_k = idxk(i);
							rowvec x_i = Covariate.row(i_k);
							rowvec w_i(2, fill::ones);
							rowvec y_i = Outcome.row(i_k);
							rowvec z_i = zz.row(i_k) - 1.0;
							w_i(1) = Treat(i_k);
							mat X(J, xcols * J, fill::zeros);
							mat W(J, nw*J, fill::zeros);
							for (int j = 0; j < J; ++j) {
								X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
								W(j, span(j*2, (j+1)*2-1)) = w_i;
							}
							vec resid = y_i.t() - arma::join_horiz(X,W) * theta - W * gamR_k - delta % z_i.t();
							loglik -= 0.5 * (vv + static_cast<double>(J)) * std::log(vv + arma::dot(resid, Sigmainv * resid));
						}
					}
					return -loglik;
				};
				double bold = std::log(nu - v0);
				double start[] = { bold };
				double xmin[] = { 0.0 };
				double ynewlo = 0.0;
				double reqmin = 1.0e-10;
				int konvge = 5;
				int kcount = 1000;
				double step[] = { 0.2 };
				int icount = 0;
				int numres = 0;
				int ifault = 0;
				nelmin(fx_vv, 1, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
				double xmax = xmin[0];
				double minll = ynewlo;

				mat cl(5,3, fill::zeros);
				vec dl(5, fill::zeros);
				double step_size = 0.1;
				vv_sampling_block:
					for (int iii=0; iii < 5; ++iii) {
						double e1 = static_cast<double>(iii-2);
						cl(iii,0) = std::pow(xmax + e1 * step_size, 2.0);
						cl(iii,1) = xmax + e1 * step_size;
						cl(iii,2) = 1.0;
						start[0] = xmax + e1 * step_size;
						dl(iii) = fx_vv(start);
					}

				for (int ni=0; ni < 5; ++ni) {
					if ((ni+1) != 3) {
						if (dl(ni) <= minll) {
							step_size *= 1.2;
							goto vv_sampling_block;
						}
					}
				}

				vec fl = arma::solve(cl.t() * cl, cl.t() * dl);
				double sigmaa = std::sqrt(0.5 / fl(0));


				double vv_prop = ::norm_rand() * sigmaa + xmax;
				start[0] = vv_prop;
				double ll_diff = fx_vv(start) - 0.5 * (std::pow(bold - xmax, 2.0) - std::pow(vv_prop - xmax, 2.0)) / std::pow(sigmaa, 2.0);
				start[0] = bold;
				ll_diff -= fx_vv(start);
				if (std::log(::unif_rand()) < ll_diff) {
					nu = v0 + std::exp(vv_prop);
					++nu_rates;
				}
				/************
				Update lambda
				************/
				double shape_lam = 0.5 * (nu + static_cast<double>(J));
				for (int k = 0; k < K; ++k) {
					uvec idxk = idxks(k);
					int n_k = idxk.n_elem;
					vec gamR_k = gamR.col(k);
					for (int i = 0; i < n_k; ++i) {
						int i_k = idxk(i);
						rowvec x_i = Covariate.row(i_k);
						rowvec w_i(2, fill::ones);
						rowvec y_i = Outcome.row(i_k);
						w_i(1) = Treat(i_k);
						mat X(J, xcols * J, fill::zeros);
						mat W(J, nw*J, fill::zeros);
						rowvec z_i = zz.row(i_k) - 1.0;

						for (int j = 0; j < J; ++j) {
							X(j, span(j*xcols, (j+1)*xcols-1)) = x_i;
							W(j, span(j*2, (j+1)*2-1)) = w_i;
						}
						vec resid = y_i.t() - arma::join_horiz(X,W) * theta - W * gamR_k - delta % z_i.t();
						lambda(i_k) = ::Rf_rgamma(shape_lam, 1.0) / (0.5 * (nu + arma::dot(resid, resid)));
					}
				}
			}

			thetaSave.col(ikeep) = theta;
			deltaSave.col(ikeep) = delta;
			lambdaSave.col(ikeep) = lambda;
			nuSave(ikeep) = nu;
			prog.increment();	
		}
	}
	

	return ListBuilder()
	.add("theta", thetaSave)
	.add("delta", deltaSave)
	.add("nu", nuSave)
	.add("lambda", lambdaSave)
	.add("nu_acceptance", nu_rates);
}