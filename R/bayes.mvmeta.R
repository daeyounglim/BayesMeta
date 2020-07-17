bayes.mvmeta <- function(Outcome, Trial, Patient, Treat, Covariate, prior=list(), mcmc=list(), verbose=FALSE) {
	if (is.matrix(Outcome)) {
		stop("Outcome must be a matrix.")
	}
	if (is.matrix(Covariate)) {
		stop("Covariate must be a matrix.")
	}
	if (length(unique(Treat)) != 2) {
		stop("Treatment should be binary.")
	}
	Treat.order <- sort(unique(Treat))
	Treat <- relabel.vec(Treat, Treat.order) - 1 # relabel the treatment numbers	

	mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
	mcvals[names(mcmc)] = mcmc
	ndiscard <- mcvals$ndiscard
	nskip <- mcvals$nskip
	nkeep <- mcvals$nkeep

	unique.trials <- unique(Trial)
	Npt <- vapply(unique.trials, function(xx) sum(Trial == xx), FUN.VALUE = integer(1))

	Trial.order <- sort(unique(Trial))
	Trial <- relabel.vec(Trial, Trial.order) # relabel the treatment numbers

	Patient.order <- sort(unique(Patient))
	Patient <- relabel.vec(Patient, Patient.order) # relabel the treatment numbers

	o <- order(Trial, Patient)
	Outcome <- Outcome[o,]
	Trial <- Trial[o]
	Patient <- Patient[o]
	Covariate <- Covariate[o,]

	priorvals <- list(a1 = 1,
					  a2 = 0.1,
					  a3 = 1,
					  a4 = 1,
					  a5 = 0.1,
					  b1 = 0.1,
					  b2 = 0.1,
					  b3 = 0.1,
					  b4 = 0.1,
					  b5 = 0.1,
					  c1 = 100,
					  c2 = 100,
					  c3 = 100,
					  d1 = )
	priorvals[names(prior)] <- prior
	nu <- priorvals$nu
	c01 <- priorvals$c01
	c02 <- priorvals$c02
}