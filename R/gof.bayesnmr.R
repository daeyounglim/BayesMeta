#' get goodness of fit 
#' @param object the output model from fitting a meta analysis/regression model
#' @param type the type of goodness of fit to compute; DIC or LPML
#' @param verbose FALSE by default; If TRUE, then progress bar will appear
#' @method gof bayesnmr
#' @export

"gof.bayesnmr" <- function(object, type="lpml", verbose=FALSE) {
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

	if (type == "dic") {
		gof <- .Call(`_BayesMeta_calc_modelfit_dic`,
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
	} else if (type == "lpml") {
		gof <- .Call(`_BayesMeta_calc_modelfit_lpml`,
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
	}

	class(gof) <- "gofnmr"
	gof
}