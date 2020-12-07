

SparseAssignTest1 <- function(sig.counts,
                              trace = 0,
                              max.mc.cores = 1) {

  set.seed(101010, kind = "L'Ecuyer-CMRG")

  sig.names <- names(sig.counts)

  some.sigs  <-
    PCAWG7::signature$genome$SBS96[ , sig.names, drop = FALSE]
  ref.genome <- attr(some.sigs, "ref.genome", exact = TRUE)
  region     <- attr(some.sigs, "region", exact = TRUE)
  if (is.null(region)) {
    message("Null region, why?")
    region <- "genome"
  }

  spect <- round(some.sigs %*% sig.counts)
  spect <-
    ICAMS::as.catalog(
      spect,
      ref.genome   = ref.genome,
      region       = region,
      catalog.type = "counts")
  nbinom.size <- 5

  m.opts <- DefaultManyOpts()
  m.opts$trace <- trace

  SA.out <- SparseAssignActivity1(spect       = spect,
                                  sigs         = some.sigs,
                                  m.opts       = m.opts,
                                  max.mc.cores = max.mc.cores
  )

  zeros <- which(SA.out < 0.5)
  if (length(zeros) > 0) {
    SA.out2 <- SA.out[-zeros]
    sig.names2 <- sig.names[-zeros]
  } else {
    SA.out2 <- SA.out
    sig.names2 <- sig.names
  }

  recon1 <- round(prop.reconstruct(some.sigs, SA.out))

  return(list(soln1       = SA.out,
              # soln2       = polish.out,
              truth       = sig.counts,
              edist1      = EDist2Spect(SA.out, sig.names, spect),
              edist1r     = EDist2SpectRounded(SA.out, sig.names, spect),
              # Was mSigBG::
              LL1         = LLHSpectrumNegBinom(spect, recon1, nbinom.size)
  ))

}


test_that("SparseAssignActivity1 (one spectrum) Test 1", {
  retval <-  SparseAssignTest1(sig.counts = c(SBS1 = 1000, SBS22 = 2000))

  testthat::expect_equal(retval$soln1,
                         c(SBS1    = 973.214846450792,
                           SBS22   = 2024.78515354921),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1),
                         14.6182783136673,
                         tolerance = 1e-2)

})


test_that("SparseAssignActivity1 (one spectrum) Test 2", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts = c(SBS3 = 10, SBS5 = 1000, SBS10a = 2000),
    trace = 0)

  testthat::expect_equal(retval$soln1,
                         c(SBS3  = 0, SBS5 = 1011.794, SBS10a = 2002.206),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1),
                         5.54163918029701,
                         tolerance = 1e-2)
})


test_that("SparseAssignActivity1 (one spectrum) Test 3", {

  # RESULTS NOT UPDATED

  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts =
      c(SBS3 = 300, SBS5 = 300, SBS4 = 300, SBS29 = 300, SBS24 = 300, SBS8 = 300),
    trace = 0)

  testthat::expect_equal(retval$soln1,
                         c(SBS3  = 299.2562,
                           SBS5  = 300.8982,
                           SBS4  = 300.0772,
                           SBS29 = 299.2562,
                           SBS24 = 299.2562,
                           SBS8  = 299.2562),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1),
                         2.845529,
                         tolerance = 1e-2)

})


test_that("SparseAssignActivity1 (one spectrum) Test 4", {

  # RESULTS NOT UPDATED

  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts =
      c(SBS3=100, SBS5=100, SBS4=100, SBS29=100, SBS24=100, SBS8=100),
    trace = 0)

  testthat::expect_equal(retval$soln1,
                         c(SBS3  = 99.83333,
                           SBS5  = 99.83333,
                           SBS4  = 99.83333,
                           SBS29 = 99.83333,
                           SBS24 = 99.83333,
                           SBS8  = 99.83333),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1), 2.921398,
                         tolerance = 1e-2)

})


test_that("SparseAssignActivity1 (one spectrum) Test 5", {

  # RESULTS NOT UPDATED

  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts = c(SBS3=30, SBS5=30, SBS4=30, SBS29=30, SBS24=30, SBS8=30),
    trace = 0)

  testthat::expect_equal(retval$soln1,
                         c(SBS3  = 0,
                           SBS5  = 63.96367,
                           SBS4  = 0,
                           SBS29 = 73.43977,
                           SBS24 = 0,
                           SBS8  = 39.59656),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1), 5.956736,
                         tolerance = 1e-2)

})
