#' Quadratic programming optimization of signature activities
#'
#' @param spectrum Mutational signature spectrum
#'      as a numeric vector or single column data frame or matrix.
#'
#' @param signatures Matrix or data frame of signatures from which
#'      reconstruct \code{spectrum}. Rows are mutation types and
#'      columns are signatures. Should have column names for
#'      interpretable results. Cannot be a vector because
#'      the column names are needed.
#'
#' @return A vector of exposures with names being the colnames from
#'   \code{signatures}.

#' Code adapted from \code{SignatureEstimation::decomposeQP}.

OptimizeExposureQP <- function(spectrum, signatures) {

  M <- spectrum / sum(spectrum)

  if (is.data.frame(signatures)) {
    P <- as.matrix(signatures)
  }
  stopifnot(is.matrix(signatures))
  P <- signatures

  # N: how many signatures should be considered
  N = ncol(P)
  if (ncol(P) == 1) {
    rr <- sum(spectrum)
    names(rr) <- colnames(P)
    stopifnot(!is.null(names(rr)))
    return(rr)
  }

  # Matrix appearing in the quadratic programming objective function
  G = t(P) %*% P

  # Constraints under which to minimize the objective function
  C <- cbind(rep(1,N), diag(N))

  # b: vector containing the values of b_0.
  b <- c(1,rep(0,N))

  # d: vector appearing in the quatric programming objective function
  d <- t(M) %*% P

  out <- quadprog::solve.QP(Dmat = G, dvec = d, Amat = C, bvec = b, meq = 1)

  # Some exposure values may be < but very close to 0;
  # Change theseto 0 and renormalize

  exposures <- out$solution
  names(exposures) <- colnames(signatures)
  exposures[exposures < 0] <- 0

  rr <- sum(spectrum) * exposures/sum(exposures)
  stopifnot(!is.null(names(rr)))
  return(rr)
}
