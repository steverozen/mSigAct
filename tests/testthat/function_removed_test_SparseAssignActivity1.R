

SparseAssignTest1 <- function(sig.counts,
                              trace = 0,
                              max.mc.cores = 1) {

  set.seed(101010, kind = "L'Ecuyer-CMRG")

  sig.names <- names(sig.counts)

  some.sigs  <-
    cosmicsig::COSMIC_v3.0$signature$GRCh37$SBS96[ , sig.names, drop = FALSE]
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

  m.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  m.opts$trace <- trace

  SA.out <- SparseAssignActivity1(spect       = spect,
                                  sigs         = some.sigs,
                                  m.opts       = m.opts,
                                  max.mc.cores = max.mc.cores
  )
  SA.out <- SA.out$proposed.assignment
  zeros <- which(SA.out < 0.5)
  if (length(zeros) > 0) {
    SA.out2 <- SA.out[-zeros]
    sig.names2 <- sig.names[-zeros]
  } else {
    SA.out2 <- SA.out
    sig.names2 <- sig.names
  }

  recon1 <- round(ReconstructSpectrum(some.sigs, SA.out, use.sig.names = TRUE))
  
  
  return(list(soln1       = SA.out[, ],
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
                         round(c(SBS1    = 973.214846450792,
                           SBS22   = 2024.78515354921)),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1),
                         14.38372,
                         tolerance = 1e-2)

})


test_that("SparseAssignActivity1 (one spectrum) Test 2", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts = c(SBS3 = 10, SBS5 = 1000, SBS10a = 2000),
    trace = 0)
  
  # RESULTS UPDATED on 9 July 2021, may need further investigation for this particular 
  # test below
  testthat::expect_equal(retval$soln1,
                         c(SBS5 = 1021, SBS10a = 1993),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1),
                         5.577088,
                         tolerance = 1e-2)
})


test_that("SparseAssignActivity1 (one spectrum) Test 3", {

  # RESULTS UPDATED

  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts =
      c(SBS3 = 300, SBS5 = 300, SBS4 = 300, SBS29 = 300, SBS24 = 300, SBS8 = 300),
    trace = 0)

  testthat::expect_equal(retval$soln1,
                         round(
                         c(SBS3  = 341.2224,
                           SBS5  = 279.0896,
                           SBS29 = 417.9223,
                           SBS24 = 350.3011,
                           SBS8  = 409.4643)),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1),
                         29.03664,
                         tolerance = 1e-2)

})


test_that("SparseAssignActivity1 (one spectrum) Test 4", {

  # RESULTS UPDATED

  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts =
      c(SBS3=100, SBS5=100, SBS4=100, SBS29=100, SBS24=100, SBS8=100),
    trace = 0)

  testthat::expect_equal(retval$soln1,
                         round(c(SBS5  = 163.8994,
                           SBS4  = 163.5813,
                           SBS24 = 153.1080,
                           SBS8  = 118.4111)),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1), 9.671682,
                         tolerance = 1e-2)

})


test_that("SparseAssignActivity1 (one spectrum) Test 5", {

  # RESULTS UPDATED

  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  retval <- SparseAssignTest1(
    sig.counts = c(SBS3=30, SBS5=30, SBS4=30, SBS29=30, SBS24=30, SBS8=30),
    trace = 0)

  testthat::expect_equal(retval$soln1,
                         round(c(SBS5  = 56.6356,
                           SBS29 = 70.1188,
                           SBS8  = 50.2455)),
                         tolerance = 1e-2)

  testthat::expect_equal(as.numeric(retval$edist1), 5.67639,
                         tolerance = 1e-2)

  unlink("Rplots.pdf")

})
