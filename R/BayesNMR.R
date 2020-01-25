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
#' @param mcmc list of MCMC-related parameters: number of burn-ins(nwarmup), number of thinning(nskip), and posterior sample size(nsave)
#' @export
bayesnmr <- function(y, sd, x, ids, iarm, groupinfo, npt, prior = list(), mcmc = list(), add.z=list(), scale.x=TRUE) {

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

	nz <- length(groupinfo)
	ns <- length(y)
	K <- length(unique(ids))
	nT <- length(unique(iarm))
	z <- matrix(0, ns, nz)
	for (j in 1:nz) {
		for (i in 1:ns) {
			if (iarm[i] %in% groupinfo[[j]]) {
				z[i,j] <- 1
			}
		}
	}
	z <- cbind(1,z)
	if (length(add.z) > 0) {
		z <- cbind(z, scale(add.z, center=TRUE, scale=TRUE))
	}
	if (scale.x) {
		x <- scale(x, center = TRUE, scale = TRUE)
	}


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
					  as.integer(nkeep))
			})

	out <- list(y = y,
				sd = sd,
				npt = npt,
				x = x,
				z = z,
				ids = ids,
				iarm = iarm,
				groupinfo = groupinfo,
				prior = priorvals,
				mcmctime = mcmctime,
				mcmc = mcvals,
				mcmc.draws = fout)
	class(out) = "bayesnmr"
	out
}
