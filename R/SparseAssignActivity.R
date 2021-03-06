#' Find known signatures that can most sparsely reconstruct each spectrum in a catalog.
#'
#' @param spectra The spectra (multiple spectra) to be reconstructed.
#'
#' @param sigs The known signatures to use in reconstruction.
#'
#' @param max.level The largest number of signatures to consider discarding
#' in the reconstruction.
#'
#' @param p.thresh The maximum p value based on which it is decided
#' to retain a signature in a reconstruction.
#'
#' @param m.opts For documentation
#'    see \code{\link{DefaultManyOpts}}.
#'
#' @param num.parallel.samples The (maximum) number of samples to run in parallel; each
#'    sample in turn can require multiple cores, as governed by
#'    \code{mc.cores.per.sample}.
#'
#' @param mc.cores.per.sample
#'    The maximum number of cores to use for each sample.
#'    On Microsoft Windows machines it is silently changed to 1.
#'
#' @return A list with the inferred exposure matrix as element \code{exposure}.
#'
#' @export

SparseAssignActivity <-
  function(spectra,
           sigs,
           max.level            = 5,
           p.thresh             = 0.05,
           m.opts               = NULL,
           num.parallel.samples = 5,
           mc.cores.per.sample  = min(20, 2^max.level)) {
    f1 <- function(i) {
      retval1 <- SparseAssignActivity1(
      spect        = spectra[ , i, drop = FALSE],
      sigs         = sigs,
      p.thresh     = p.thresh,
      m.opts       = m.opts,
      max.level    = max.level,
      max.mc.cores = mc.cores.per.sample)

    return(retval1)
  }

  if (is.null(m.opts)) m.opts <- DefaultManyOpts()

  num.parallel.samples <- Adj.mc.cores(num.parallel.samples)

  retval <- parallel::mclapply(1:ncol(spectra),
                               f1,
                               mc.cores = num.parallel.samples)
  check.mclapply.result(
    retval, "SparseAssignActivity", colnames(spectra))

  other.info <- lapply(retval, attributes)
  retval2 <- matrix(unlist(retval), ncol = length(retval))
  colnames(retval2) <- colnames(spectra)
  rownames(retval2) <- colnames(sigs)

  return(list(exposure = retval2, other.info = other.info))
}

