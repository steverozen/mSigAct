#' Test whether a given signature is plausibly present in a spectrum.
#'
#' For backward compatibility. See also \code{\link{AnySigSubsetPresent}}.
#'
#' @param spectrum The spectrum to analyze.
#'
#' @param sigs A catalog of signatures from which to choose.
#'
#' @param target.sig.index The index of the signature the presence
#' of which we want to test. It can also be the signature id (e.g. "SBS22").
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'    
#' @param seed Random seed; set this to get reproducible results. (The
#'   numerical optimization is in two phases; the first, global phase
#'   might rarely find different optima depending on the random
#'   seed.)
#' 
#' @return A list of elements:
#' * loglh.with: The maximum log likelihood of the reconstructed spectrum using
#' all the signatures.
#'
#' * loglh.without: The maximum log likelihood of the reconstructed spectrum
#' without the target signature.
#'
#' * statistic: Likelihood ratio test statistic.
#'
#' * chisq.p: P-value of the likelihood ratio test. The null hypothesis is we
#' can plausibly reconstruct the spectrum without the target signature.
#'
#' * exp.with: The exposure using all the signatures which generates the maximum
#' log likelihood \code{loglh.with}.
#'
#' * exp.without: The exposure not using the target signature which generates
#' the maximum log likelihood \code{loglh.without}.
#' 
#' @md
#' 
#' @keywords internal

SignaturePresenceTest1 <- function(
  spectrum, sigs, target.sig.index, m.opts = DefaultManyOpts(), seed = NULL) {
  
  if (!is.numeric(target.sig.index)) {
    if (!target.sig.index %in% colnames(sigs)) {
      stop("The signature id specified by target.sig.index ", target.sig.index,
           " is not found in column names of sigs")
    }
    target.sig.index <- which(colnames(sigs) == target.sig.index)
  }
  
  if (!is.null(seed)) set.seed(seed, kind = "L'Ecuyer-CMRG")
  
  ret.with <- OptimizeExposure(spectrum  = spectrum,
                               sigs   = sigs,
                               m.opts = m.opts)
  loglh.with <- ret.with$loglh

  ret.without <- OptimizeExposure(spectrum = spectrum,
                                  sigs   = sigs[ ,-target.sig.index],
                                  m.opts = m.opts)
  loglh.without <- ret.without$loglh

  statistic <- 2 * (loglh.with - loglh.without)
  chisq.p <- stats::pchisq(statistic, 1, lower.tail = FALSE)
  if (m.opts$trace > 0)
    message("statistic  = ", statistic, "\nchisq p = ", chisq.p)

  list(loglh.with              = loglh.with,
       loglh.without           = loglh.without,
       statistic               = statistic,
       chisq.p                 = chisq.p,
       exp.with                = ret.with$exposure,
       exp.without             = ret.without$exposure)
       #everything.else.with    = ret.with$everything.else,
       #everything.else.without = ret.without$everything.else)
}

#' Framework for testing \code{\link{SignaturePresenceTest1}}.
#'
#' @keywords internal
#'
TestSignaturePresenceTest1 <-
  function(sig.counts,
           input.sigs = PCAWG7::signature$genome$SBS96,
           trace      = 0,
           m.opts     = NULL) {

    if (!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
                          quietly = TRUE)) {
      stop("Please install Bioconductor library ",
           "BSgenome.Hsapiens.1000genomes.hs37d5")
    }

    sig.names <- names(sig.counts)

  if (sum(sig.names %in% colnames(input.sigs)) == 0) {
    stop("sig.names not all in input.sigs")
  }

  some.sigs  <- input.sigs[ , sig.names, drop = FALSE]

  ref.genome <- attr(some.sigs, "ref.genome", exact = TRUE)
  region     <- attr(some.sigs, "region", exact = TRUE)
  if (is.null(region)) {
    message("Null region, why?")
    region <- "genome"
  }

  spectrum <- round(some.sigs %*% sig.counts)
  spectrum <-
    ICAMS::as.catalog(
      object       = spectrum,
      ref.genome   = ref.genome,
      region       = region,
      catalog.type = "counts")

  if (is.null(m.opts)) {
    m.opts <- DefaultManyOpts()
  }
  m.opts$trace <- trace

  retval <- SignaturePresenceTest1(
    spectrum         = spectrum,
    sigs             = some.sigs,
    target.sig.index = 1,
    m.opts           = m.opts)

  return(retval)
}


# https://cran.r-project.org/web/packages/DirichletReg/DirichletReg.pdf


