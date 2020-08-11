#' Fit Bayesian multivariate skew meta-regression models for individual patient data
#' 
#' This is a function for running the Markov chain Monte Carlo algorithm for the BFSHMVMR Model. The first six arguments are required.
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param Outcome the responses for each patient in every trial
#' @param Trial the identifier for a trial; it will internally be relabeled to contiguous/consecutive natural numbers starting from 1
#' @param Patient the patient number in each trial
#' @param Treat treatment applied; it will internally be relabeled to 0 and 1
#' @param Covariate covariates for fixed effects
#' @param tau hyperparameter which determines the level of skewness
#' @param prior list of hyperparameters; when not given, algorithm will run in default setting
#' @param mcmc list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @return a dataframe with input arguments, posterior samples, Metropolis algorithm acceptance rates, etc
#' @examples
#' \dontrun{
#' data(df)
#' groupInfo <- list(c(0,1), c(2,3), c(4)) # define the variance structure
#' x <- df[,6:10]
#' fit <- bayes.mvmeta(df$y, df$sd, x, df$ids, df$iarm, df$npt, groupInfo,
#' 			prior = list(c01=1.0e05, c02=4, nu=3),
#' 			mcmc=list(ndiscard=2500,nskip=1,nkeep=10000))
#' }
#' @export
bayes.mvmeta <- function(Outcome, Trial, Patient, Treat, Covariate, tau, prior=list(), mcmc=list(), verbose=FALSE) {
	if (!is.matrix(Outcome)) {
		stop("Outcome must be a matrix.")
	}
	if (!is.matrix(Covariate)) {
		stop("Covariate must be a matrix.")
	}
	unique.treat <- unique(Treat)
	if (length(unique.treat) != 2) {
		stop("Treatment should be binary.")
	}
	Treat.order <- sort(unique.treat)
	Treat <- relabel.vec(Treat, Treat.order) - 1 # relabel the treatment numbers	

	mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
	mcvals[names(mcmc)] = mcmc
	ndiscard <- mcvals$ndiscard
	nskip <- mcvals$nskip
	nkeep <- mcvals$nkeep

	unique.trials <- unique(Trial)
	K <- length(unique.trials)
	J <- ncol(Outcome)
	# Npt <- vapply(unique.trials, function(xx) sum(Trial == xx), FUN.VALUE = integer(1))

	Trial.order <- sort(unique.trials)
	Trial <- relabel.vec(Trial, Trial.order) # relabel the treatment numbers

	Patient.order <- sort(unique(Patient))
	Patient <- relabel.vec(Patient, Patient.order) # relabel the treatment numbers

	o <- order(Trial, Patient)
	Outcome <- Outcome[o,]
	Trial <- Trial[o]
	Patient <- Patient[o]
	Covariate <- Covariate[o,]



	priorvals <- list(a0 = 1,
					  b0 = 0.1,
					  v0 = 3,
					  c1 = 100,
					  c2 = 100,
					  c3 = 100,
					  d0 = J + 0.1,
					  S0 = diag(0.1, J),
					  d1 = 2.1,
					  S1 = diag(0.1, 2))
	priorvals[names(prior)] <- prior
	a0 <- priorvals$a0
	b0 <- priorvals$b0
	v0 <- priorvals$v0
	c1 <- priorvals$c1
	c2 <- priorvals$c2
	c3 <- priorvals$c3
	d0 <- priorvals$d0
	S0 <- priorvals$S0
	d1 <- priorvals$d1
	S1 <- priorvals$S1
	

	mcmctime <- system.time({		
		fit <- .Call(`_BayesMeta_BMVSMR`,
					as.matrix(Outcome),
					as.matrix(Covariate),
					as.double(Treat),
					as.integer(Patient),
					as.integer(Trial),
					as.double(a0),
					as.double(b0),
					as.double(c1),
					as.double(c2),
					as.double(c3),
					as.double(d0),
					as.matrix(S0),
					as.double(d1),
					as.matrix(S1),
					as.double(v0),
					as.double(tau),
					as.integer(K),
					as.integer(ndiscard),
					as.integer(nskip),
					as.integer(nkeep),
					as.logical(verbose))
	})
	fit
}