#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param y aggregate mean of the responses for each arm of each study
#' @param sd standard deviation of the responses for each arm of each study
#' @param x aggregate covariates for the mean component
#' @param ids study number in integers
#' @param iarm arm number of each study
#' @param groupinfo list of grouping information; the control(baseline) group must start from 0; the aggregate covariates 'z' explaining the variance of the random effect of the t-th treatment will be construct based on this grouping information
#' @param npt number of observations per trial
#' @param nT number of treatments
#' @param prior list of hyperparameters; when not given, algorithm will run in default setting
#' @param mcmc list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param add.z additional covariates other than the grouping vectors that should be column-concatenated to 'z'. This should have the same number of rows as 'y', and 'x'
#' @param scale.x logical variable for scaling x. Default to TRUE. If not scaled, the gamma[1] (different than gam in the function) cannot be interpreted as placebo
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @param init initial values for beta (ns + nT dimensional) and phi. Dimensions must be conformant.
#' @export
bayesnmr <- function(y, sd, x, ids, iarm, npt, groupinfo=list(), prior = list(), mcmc = list(), add.z=list(), scale.x=TRUE, verbose=FALSE, init=list()) {

	mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
	mcvals[names(mcmc)] = mcmc
	ndiscard <- mcvals$ndiscard
	nskip <- mcvals$nskip
	nkeep <- mcvals$nkeep


	priorvals <- list(nu = 20, c01 = 1.0e05, c02 = 4)
	priorvals[names(prior)] <- prior
	nu <- priorvals$nu
	c01 <- priorvals$c01
	c02 <- priorvals$c02

	if (min(unique(iarm)) != 0) {
		stop("The treatments should start from 0. Please adjust accordingly.")
	}
	nx <- ncol(x)
	nz <- length(groupinfo)
	ns <- length(y)
	K <- length(unique(ids))
	nT <- length(unique(iarm))
	z <- matrix(0, ns, nz)
	if (nz > 0) {
		for (j in 1:nz) {
			for (i in 1:ns) {
				if (iarm[i] %in% groupinfo[[j]]) {
					z[i,j] <- 1
				}
			}
		}
	}
	# z <- cbind(1,z)
	if (length(add.z) > 0) {
		z <- cbind(z, scale(add.z, center=TRUE, scale=TRUE))
	}
	if (scale.x) {
		x <- scale(x, center = TRUE, scale = TRUE)
	}
	init_final <- list(beta = numeric(nx+nT), phi = numeric(ncol(z)), sig2 = rep(1, ns))
	init_final[names(init)] <- init


	mcmctime <- system.time({
				fout <- .Call(`_BayesMeta_BayesNMR`,
					  as.double(y),
					  as.double(sd),
					  as.matrix(x),
					  as.matrix(z),
					  as.integer(ids),
					  as.integer(iarm),
					  as.double(npt),
					  as.double(nu),
					  as.double(1/c01),
					  as.double(1/c02),
					  as.integer(K),
					  as.integer(nT),
					  as.integer(ndiscard),
					  as.integer(nskip),
					  as.integer(nkeep),
					  as.logical(verbose),
					  as.double(init_final$beta),
					  as.double(init_final$phi),
					  as.double(init_final$sig2))
			})

	out <- list(y = y,
				sd = sd,
				npt = npt,
				x = x,
				z = z,
				ids = ids,
				iarm = iarm,
				K = K,
				nT = nT,
				groupinfo = groupinfo,
				prior = priorvals,
				mcmctime = mcmctime,
				mcmc = mcvals,
				mcmc.draws = fout)
	class(out) <- "bayesnmr"
	out
}
