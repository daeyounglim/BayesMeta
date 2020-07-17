#ifndef BAYESMETA_RANDOM_H
#define BAYESMETA_RANDOM_H
#include <stdio.h>
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include <Rdefines.h>

using namespace Rcpp;
using namespace R;
using namespace arma;

arma::mat rtmvn_gibbs(const int& n, const int& p, const arma::vec& Mean, const arma::mat& Sigma_chol,
                      const arma::mat& R, const arma::vec& a, const arma::vec& b, const arma::vec& z);

arma::mat rtmvnCpp(const int& n, const arma::vec& Mean, const arma::mat& Sigma, const arma::mat& D,
                   const arma::vec& lower, const arma::vec& upper, const arma::vec& init);