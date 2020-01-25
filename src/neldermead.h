#ifndef BAYESMETA_NELDERMEAD_H
#define BAYESMETA_NELDERMEAD_H
# include <cmath>
# include <RcppArmadillo.h>
#include <Rdefines.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;

void nelmin ( double fn ( const arma::vec& x ),
              int n,
              arma::vec& start,
              arma::vec& xmin, 
              double& ynewlo,
              double& reqmin,
              arma::vec& step,
              int konvge,
              int kcount, 
              int& icount,
              int& numres,
              int& ifault );

#endif
