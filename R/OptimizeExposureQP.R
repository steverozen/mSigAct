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
#' @export
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

#' Bootstrap \code{\link{OptimizeExposureQP}} and filter exposures by confidence intervals
#'
#' @export
#'
#' @inheritParams OptimizeExposureQP
#'
#' @param num.replicates Number of bootstrap replicates.
#'
#' @param conf.interval Discard signatures with \code{conf.int} that overlaps 0.
#'
#' @param mc.cores The maximum number of cores to use.
#'   On MS Windows machines it defaults to 1.

OptimizeExposureQPBootstrap <- function(spectrum,
                                        signatures,
                                        num.replicates = 10000,
                                        conf.int       = 0.95,
                                        mc.cores       = 10) {

  mc.cores <- Adj.mc.cores(mc.cores) # Set to 1 if OS is MS Windows
  ss <- sum(spectrum)
  spectrum.as.probs <- spectrum/ss

  s2 <- stats::rmultinom(
    n = num.replicates, size = ss, prob = spectrum.as.probs)

  my.fn <- function(ii) {
    ex <- OptimizeExposureQP(spectrum = s2[ , ii], signatures)
    return(ex)
  }

  rr <- parallel::mclapply(
    X = 1:num.replicates, FUN = my.fn, mc.cores = mc.cores)
  rr2 <- do.call(cbind, rr)
  ex.lower.ci <- apply(rr2, MARGIN = 1, FUN = stats::quantile, probs =  conf.int / 2)

  while (any(ex.lower.ci == 0)) {
    signatures <- signatures[ , which(ex.lower.ci > 0), drop = FALSE]
    rr <- parallel::mclapply(
      X = 1:num.replicates, FUN = my.fn, mc.cores = mc.cores)
    rr2 <- do.call(cbind, rr)
    ex.lower.ci <- apply(rr2, MARGIN = 1, FUN = stats::quantile, probs =  conf.int / 2)
  }

  exp <- OptimizeExposureQP(spectrum, signatures)
  return(exp)
}

if (FALSE) {
  foo <- OptimizeExposureQPBootstrap(spectrum = PCAWG7::spectra$PCAWG$SBS96[ , 1],
                              signatures = PCAWG7::signature$genome$SBS96,
                              mc.cores = 50,
                              num.replicates = 10000)

}


