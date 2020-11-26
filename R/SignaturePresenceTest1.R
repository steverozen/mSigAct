#' Test whether a given signature is plausibly present in a spectrum.
#'
#' Use \code{\link{AnySigSubsetPresent}}instead.
#'
#' @param spectrum The spectrum to analyze
#'
#' @param sigs A catalog of signatures from which to choose
#'
#' @param target.sig.index The index of the signature the presence
#' of which we want to test.
#'
#' @param eval_f See \code{\link[nloptr]{nloptr}}.
#'
#' @param m.opts For documentation
#'    see \code{\link{DefaultManyOpts}}.
#'
#' @keywords internal

SignaturePresenceTest1 <- function(
  spectrum, sigs, target.sig.index, m.opts, eval_f) {

  ret.with <- OptimizeExposure(spectrum  = spectrum,
                               sigs   = sigs,
                               m.opts = m.opts,
                               eval_f = eval_f)
  loglh.with <- ret.with$loglh

  ret.without <- OptimizeExposure(spectrum = spectrum,
                                  sigs   = sigs[ ,-target.sig.index],
                                  eval_f = eval_f,
                                  m.opts = m.opts)
  loglh.without <- ret.without$loglh

  statistic <- 2 * (loglh.with - loglh.without)
  chisq.p <- stats::pchisq(statistic, 1, lower.tail = FALSE)
  if (m.opts$trace > 0)
    message("statistic  = ", statistic, "\nchisq p = ", chisq.p)

  list(with                    = loglh.with,
       without                 = loglh.without,
       statistic               = statistic,
       chisq.p                 = chisq.p,
       exp.with                = ret.with$exp,
       exp.without             = ret.without$exp,
       everything.else.with    = ret.with$everything.else,
       everything.else.without = ret.without$everything.else)
}

#' Framework for testing \code{\link{SignaturePresenceTest1}}.
#'
#' @keywords internal
#'
TestSignaturePresenceTest1 <-
  function(sig.counts,
           input.sigs = PCAWG7::signature$genome$SBS96,
           trace      = 0,
           eval_f     = ObjFnBinomMaxLHRound,
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
    eval_f           = eval_f,
    m.opts           = m.opts)

  return(retval)
}


# https://cran.r-project.org/web/packages/DirichletReg/DirichletReg.pdf


