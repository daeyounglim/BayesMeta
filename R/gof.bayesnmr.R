"gof.bayesnmr" <- function(object, verbose=FALSE) {
	y <- object$y
	npt <- object$npt
	x <- object$x
	z <- object$z
	ids <- object$ids
	iarm <- object$iarm
	K <- object$K
	nT <- object$nT
	nkeep <- object$mcmc$nkeep
	nu <- object$prior$nu

	gof <- .Call(`_BayesMeta_calc_modelfit`,
				 as.double(y),
				 as.matrix(x),
				 as.matrix(z),
				 as.integer(ids),
				 as.integer(iarm),
				 as.double(npt),
				 as.double(nu),
				 as.matrix(object$mcmc.draws$beta),
				 as.matrix(object$mcmc.draws$sig2),
				 as.matrix(object$mcmc.draws$phi),
				 as.matrix(object$mcmc.draws$lam),
				 as.array(object$mcmc.draws$Rho),
				 as.integer(K),
				 as.integer(nT),
				 as.integer(nkeep),
				 as.logical(verbose))
	class(gof) <- "gofnmr"
	gof
}