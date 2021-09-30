#' Bootstrap \code{\link{OptimizeExposureQP}} and filter exposures by confidence intervals
#'
#' @export
#'
#' @inheritParams OptimizeExposureQP
#'
#' @param num.replicates Number of bootstrap replicates.
#'
#' @param conf.int Discard signatures with \code{conf.int} that overlaps 0.
#'
#' @param mc.cores The maximum number of cores to use.
#'   On MS Windows machines it defaults to 1.
#'   
#' @param seed Random seed; set this to get reproducible results.
#'
#'   #' @return
#' A list with elements \describe{
#' \item{\code{exposure}}{The vector of exposures that generated \code{loglh},
#'    i.e. the number of mutations ascribed to each signature. The names of
#'    \code{exposure} are a subset of the \code{colnames(signatures)}.}
#' \item{\code{euclidean.dist}}{The final value of the objective function.}
#' \item{\code{cosine.sim}}{The cosine similarity between \code{spectrum} and
#'              the reconstruction based on \code{signatures}.}
#' }
#' 
#' If the spectrum has 0 mutations, no bootstrapping is done, and in the
#' return value
#' all \code{signaures} have 0 exposures, \code{euclidian.dist} is 0,
#' and \code{cosine.sim} is \code{NaN}.
#'

OptimizeExposureQPBootstrap <- function(spectrum,
                                        signatures,
                                        num.replicates = 10000,
                                        conf.int       = 0.95,
                                        mc.cores       = 10,
                                        seed           = NULL) {
  
  ss <- sum(spectrum)
  if (ss == 0) {
    return(list(exposure = rep(0, ncol(signatures)), 
                euclidean.dist = 0, 
                cosine.sim = NaN))
  }
  if (!is.null(seed)) set.seed(seed, kind = "L'Ecuyer-CMRG")
  mc.cores <- Adj.mc.cores(mc.cores) # Set to 1 if OS is MS Windows
  
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
  
  recon.spect <- as.vector(ReconstructSpectrum(sigs = signatures, exp = exp))
  if (length(recon.spect) != length(spectrum)) {
    msg <- paste("OptimizeExposureQPBootstrap:\n",
                 "length(recon.spect) =", length(recon.spect),
                 "length(spectrum) =", length(spectrum))
    stop(msg)
    
  }
  if (is.matrix(spectrum)) spectrum <- spectrum[, 1]
  rb <- rbind(spectrum, recon.spect)
  edist <- philentropy::distance(rb, method = "euclidean", test.na = FALSE)
  csim  <- philentropy::distance(rb, method = "cosine", test.na = FALSE)
  
  return(list(exposure = exp, euclidean.dist = edist, cosine.sim = csim))
}



test_that("OptimizeExposureQPBootstrap", {
  set.seed(999, kind = "L'Ecuyer-CMRG")
  rr <- OptimizeExposureQPBootstrap(spectrum = PCAWG7::spectra$PCAWG$SBS96[ , 1],
                                    signatures = PCAWG7::COSMIC.v3.0$signature$genome$SBS96[ , 1:5],
                                    mc.cores = 1,
                                    num.replicates = 100)
  expect_equal(rr$exposure,
               c(SBS1 = 1384.08972804323, SBS2 = 1684.34516525904,
                 SBS3 = 6960.32503149651,  SBS4 = 1196.35199280255,
                 SBS5 = 3674.88808239867) )
  expect_equal(rr$cosine.sim, c(cosine = 0.933303527760909))
  expect_equal(rr$euclidean.dist, c(euclidean = 760.830508126611))
})
