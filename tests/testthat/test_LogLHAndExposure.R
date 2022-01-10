PrepOneSynSpectrum <- function(sig.counts,
                               # input.sigs = PCAWG7::signature$genome$SBS96
                               input.sigs = PCAWG7::COSMIC.v3.0$signature$genome$SBS96
                               ) {

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

  return(list(sigs = some.sigs, spec = spectrum))
}

TestOneLLHetc <- function(sig.counts,
                          # input.sigs   = PCAWG7::signature$genome$SBS96,
                          input.sigs = cosmicsig::COSMIC_v3.0$signature$GRCh37$SBS96,
                          trace = 0) {

  test.data <- PrepOneSynSpectrum(sig.counts = sig.counts,
                                  input.sigs = input.sigs)


  m.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  m.opts$trace <- trace

  retval <- OptimizeExposure(
    spect       = test.data$spec,
    sigs        = test.data$sigs,
    m.opts      = m.opts)

  new.rec <-
    mSigAct:::prop.reconstruct(exp = retval$exposure, sigs = test.data$sigs)

  xx <- rbind(test.data$spec[ ,1], round(new.rec)[, 1])
  edist <- stats::dist(xx, method = "euclidean")

  return(c(retval, list(m.opts = m.opts, edist = edist)))
}

test_that("test-LogLHAndExposure.R LogLHAndExposure 1", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  set.seed(101010, kind = "L'Ecuyer-CMRG")
  retval <- TestOneLLHetc(sig.counts = input)
  testthat::expect_equal(retval$loglh, -225.1579, tolerance = 1e-5)
  testthat::expect_equal(retval$exposure,
                         c(SBS1 = 974.8783, SBS22 = 2024.1217))
  testthat::expect_equal(retval$m.opts$nbinom.size, 5)
  testthat::expect_equal(retval$m.opts$global.opts$algorithm,
                         "NLOPT_GN_DIRECT")
})
