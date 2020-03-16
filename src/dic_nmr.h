#ifndef BAYESMETA_DICNMR_H
#define BAYESMETA_DICNMR_H
#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <iterator>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rdefines.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
// [[Rcpp::depends(RcppArmadillo, RcppProgress, BH)]]
using namespace arma;

/********************************
negative log-likelihood of lambda
for Nelder-Mead method
********************************/
class NegLogLikLam
{

private: 
	const double& nu;
	const arma::vec& resid_k;
	const arma::mat& ZEREZ_S;

public:
	arma::vec start;
	arma::vec xmin;
	double ynewlo;
	double reqmin;
	arma::vec step;
	int konvge;
	int kcount;
	int icount;
	int numres;
	int ifault;

	NegLogLikLam(const double& nu_,
				 const arma::vec& resid_k_,
				 const arma::mat& ZEREZ_S_) : nu(nu_), resid_k(resid_k_), ZEREZ_S(ZEREZ_S_) {
		start = arma::zeros<arma::vec>(1);
		xmin = arma::zeros<arma::vec>(1);
		xmin(0) = 0.0;
		ynewlo = 0.0;
		reqmin = 1.0e-20;
		step = arma::zeros<arma::vec>(1);
		step(0) = 0.2;
		konvge = 5;
		kcount = 1000;
		icount = 0;
		numres = 0;
		ifault = 0;
	}
	
	double fn(const arma::vec& lam) const {
		double lam_k = lam(0);
		double eta_k = std::log(lam_k);
		double loglik = 0.5 * nu * (eta_k - lam_k) - eta_k;
		double logdet_val;
		double logdet_sign;
		arma::log_det(logdet_val, logdet_sign, ZEREZ_S);
		loglik += -0.5 * logdet_val - 0.5 * arma::accu(resid_k % arma::solve(ZEREZ_S, resid_k));
		return -loglik;
	}

	void optimize(const double& start_) {
		using namespace arma;
		int n = 1;

		double ccoeff = 0.5;
		double del;
		double dn;
		double dnn;
		double ecoeff = 2.0;
		double eps = 0.001;
		int i;
		int ihi;
		int ilo;
		int j;
		int jcount;
		int l;
		int nn;

		double rcoeff = 1.0;
		double rq;
		double x;
		double y2star;
		double ylo;
		double ystar;
		double z;
		//
		//  Check the input parameters.
		//
		if ( reqmin <= 0.0 )
		{
			Rf_error("NELMIN ERROR: reqmin, n, or konvge has an illegal value.");
		}

		if ( n < 1 )
		{
			Rf_error("NELMIN ERROR: reqmin, n, or konvge has an illegal value.");
		}

		if ( konvge < 1 )
		{
			Rf_error("NELMIN ERROR: reqmin, n, or konvge has an illegal value.");
		}

		vec p(n*(n+1));
		vec p2star(n);
		vec pbar(n);
		vec pstar(n);
		vec y(n+1);

		icount = 0;
		numres = 0;

		jcount = konvge; 
		dn = static_cast<double>( n );
		nn = n + 1;
		dnn = static_cast<double>( nn );
		del = 1.0;
		rq = reqmin * dn;
		//
		//  Initial or restarted loop.
		//
		for ( ; ; )
		{
		for ( i = 0; i < n; ++i )
		{ 
		  p(i+n*n) = start(i);
		}
		y(n) = fn ( start );
		++icount;

		for ( j = 0; j < n; ++j )
		{
		  x = start(j);
		  start(j) += step(j) * del;
		  for ( i = 0; i < n; ++i )
		  {
		    p(i+j*n) = start(i);
		  }
		  y(j) = fn ( start );
		  ++icount;
		  start(j) = x;
		}
		//                    
		//  The simplex construction is complete.
		//                    
		//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
		//  the vertex of the simplex to be replaced.
		//                
		ylo = y(0);
		ilo = 0;

		for ( i = 1; i < nn; ++i )
		{
		  if ( y(i) < ylo )
		  {
		    ylo = y(i);
		    ilo = i;
		  }
		}
		//
		//  Inner loop.
		//
		for ( ; ; )
		{
		  if ( kcount <= icount )
		  {
		    break;
		  }
		  ynewlo = y(0);
		  ihi = 0;

		  for ( i = 1; i < nn; ++i )
		  {
		    if ( ynewlo < y(i) )
		    {
		      ynewlo = y(i);
		      ihi = i;
		    }
		  }
		//
		//  Calculate PBAR, the centroid of the simplex vertices
		//  excepting the vertex with Y value YNEWLO.
		//
		  for ( i = 0; i < n; ++i )
		  {
		    z = 0.0;
		    for ( j = 0; j < nn; ++j )
		    { 
		      z = z + p(i+j*n);
		    }
		    z = z - p(i+ihi*n);  
		    pbar(i) = z / dn;
		  }
		//
		//  Reflection through the centroid.
		//
		  for ( i = 0; i < n; ++i )
		  {
		    pstar(i) = pbar(i) + rcoeff * ( pbar(i) - p(i+ihi*n) );
		  }
		  ystar = fn ( pstar );
		  ++icount;
		//
		//  Successful reflection, so extension.
		//
		  if ( ystar < ylo )
		  {
		    for ( i = 0; i < n; ++i )
		    {
		      p2star(i) = pbar(i) + ecoeff * ( pstar(i) - pbar(i) );
		    }
		    y2star = fn ( p2star );
		    ++icount;
		//
		//  Check extension.
		//
		    if ( ystar < y2star )
		    {
		      for ( i = 0; i < n; ++i )
		      {
		        p(i+ihi*n) = pstar(i);
		      }
		      y(ihi) = ystar;
		    }
		//
		//  Retain extension or contraction.
		//
		    else
		    {
		      for ( i = 0; i < n; ++i )
		      {
		        p(i+ihi*n) = p2star(i);
		      }
		      y(ihi) = y2star;
		    }
		  }
		//
		//  No extension.
		//
		  else
		  {
		    l = 0;
		    for ( i = 0; i < nn; ++i )
		    {
		      if ( ystar < y(i) )
		      {
		        l = l + 1;
		      }
		    }

		    if ( 1 < l )
		    {
		      for ( i = 0; i < n; ++i )
		      {
		        p(i+ihi*n) = pstar(i);
		      }
		      y(ihi) = ystar;
		    }
		//
		//  Contraction on the Y(IHI) side of the centroid.
		//
		    else if ( l == 0 )
		    {
		      for ( i = 0; i < n; ++i )
		      {
		        p2star(i) = pbar(i) + ccoeff * ( p(i+ihi*n) - pbar(i) );
		      }
		      y2star = fn ( p2star );
		      ++icount;
		//
		//  Contract the whole simplex.
		//
		      if ( y(ihi) < y2star )
		      {
		        for ( j = 0; j < nn; ++j )
		        {
		          for ( i = 0; i < n; ++i )
		          {
		            p(i+j*n) = ( p(i+j*n) + p(i+ilo*n) ) * 0.5;
		            xmin(i) = p(i+j*n);
		          }
		          y(j) = fn ( xmin );
		          ++icount;
		        }
		        ylo = y(0);
		        ilo = 0;

		        for ( i = 1; i < nn; ++i )
		        {
		          if ( y(i) < ylo )
		          {
		            ylo = y(i);
		            ilo = i;
		          }
		        }
		        continue;
		      }
		//
		//  Retain contraction.
		//
		      else
		      {
		        for ( i = 0; i < n; ++i )
		        {
		          p(i+ihi*n) = p2star(i);
		        }
		        y(ihi) = y2star;
		      }
		    }
		//
		//  Contraction on the reflection side of the centroid.
		//
		    else if ( l == 1 )
		    {
		      for ( i = 0; i < n; ++i )
		      {
		        p2star(i) = pbar(i) + ccoeff * ( pstar(i) - pbar(i) );
		      }
		      y2star = fn ( p2star );
		      ++icount;
		//
		//  Retain reflection?
		//
		      if ( y2star <= ystar )
		      {
		        for ( i = 0; i < n; ++i )
		        {
		          p(i+ihi*n) = p2star(i);
		        }
		        y(ihi) = y2star;
		      }
		      else
		      {
		        for ( i = 0; i < n; ++i )
		        {
		          p(i+ihi*n) = pstar(i);
		        }
		        y(ihi) = ystar;
		      }
		    }
		  }
		//
		//  Check if YLO improved.
		//
		  if ( y(ihi) < ylo )
		  {
		    ylo = y(ihi);
		    ilo = ihi;
		  }
		  --jcount;

		  if ( 0 < jcount )
		  {
		    continue;
		  }
		//
		//  Check to see if minimum reached.
		//
		  if ( icount <= kcount )
		  {
		    jcount = konvge;

		    z = 0.0;
		    for ( i = 0; i < nn; ++i )
		    {
		      z = z + y(i);
		    }
		    x = z / dnn;

		    z = 0.0;
		    for ( i = 0; i < nn; ++i )
		    {
		      z = z + pow ( y(i) - x, 2 );
		    }

		    if ( z <= rq )
		    {
		      break;
		    }
		  }
		}
		//
		//  Factorial tests to check that YNEWLO is a local minimum.
		//
		for ( i = 0; i < n; ++i )
		{
		  xmin(i) = p(i+ilo*n);
		}
		ynewlo = y(ilo);

		if ( kcount < icount )
		{
		  ifault = 2;
		  break;
		  // Rf_error("NELMIN ERROR: iteration terminated because KCOUNT was exceeded without convergence.");
		}

		ifault = 0;

		for ( i = 0; i < n; ++i )
		{
		  del = step(i) * eps;
		  xmin(i) += del;
		  z = fn ( xmin );
		  ++icount;
		  if ( z < ynewlo )
		  {
		    ifault = 2;
		    break;
		    // Rf_error("NELMIN ERROR: iteration terminated because KCOUNT was exceeded without convergence.");
		  }
		  xmin(i) += - del - del;
		  z = fn ( xmin );
		  ++icount;
		  if ( z < ynewlo )
		  {
		    ifault = 2;
		    break;
		    // Rf_error("NELMIN ERROR: iteration terminated because KCOUNT was exceeded without convergence.");
		  }
		  xmin(i) += del;
		}

		if ( ifault == 0 )
		{
		  break;
		}
		//
		//  Restart the procedure.
		//
		for ( i = 0; i < n; ++i )
		{
		  start(i) = xmin(i);
		}
		del = eps;
		++numres;
		}

		return;
	}
};

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
						 const arma::mat& lams,
						 const arma::cube& Rhos,
						 const int& K,
						 const int& nT,
						 const int& nkeep,
						 const bool verbose);

#endif
