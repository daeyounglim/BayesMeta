#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param A matrix
#' @param v vector to add
#' @export
testfun <- function(A, v) {
	mcmctime <- system.time({
				fout <- .Call(`_BayesMeta_testfun`,
					  as.matrix(A),
					  as.vector(v))
			})
	return(fout)
}
