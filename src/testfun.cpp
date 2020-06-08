#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List testfun(const arma::mat& A, const arma::vec& v) {
	arma::mat B = A;
	B.diag() += v;
	arma::mat C = A;
	C += arma::diagmat(v);
	return Rcpp::List::create(Rcpp::Named("inplace") = B, Rcpp::Named("diagmat") = C);
}

