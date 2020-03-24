context("test-SparseAssignActivity1.R")

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
                                  eval_f       = ObjFnBinomMaxLHMustRound,
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
  
  if (FALSE) {
    
    polish.out <- Polish(exp       = SA.out2,
                         sig.names = sig.names2,
                         spect     = spect)
  }
  
  recon1 <- round(prop.reconstruct(some.sigs, SA.out))
  
  if (FALSE) {
    recon2 <-
      round(
        prop.reconstruct(
          PCAWG7::signature$genome$SBS96[ , sig.names2, drop = FALSE],
          polish.out))
  }
  
  return(list(soln1       = SA.out,
              # soln2       = polish.out,
              truth       = sig.counts,
              edist1      = EDist2Spect(SA.out, sig.names, spect),
              edist1r     = EDist2SpectRounded(SA.out, sig.names, spect),
              LL1         = LLHSpectrumNegBinom(spect, recon1, nbinom.size) #,
              # LL2         = LLHSpectrumNegBinom(spect, recon2, nbinom.size),
              # edist2      = EDist2Spect(polish.out, sig.names2, spect),
              # ed8st2r     = EDist2SpectRounded(polish.out, sig.names2, spect)
              #, input.spect = spect
  ))
  
}


test_that("SparseAssignActivity1 (one spectrum) Test 1", {
  retval <-  SparseAssignTest1(sig.counts = c(SBS1 = 1000, SBS22 = 2000))
  
  testthat::expect_equal(retval$soln1,
                         c(SBS1    = 999.2467, 
                           SBS22   = 1998.7533),
                         tolerance = 1e-2)
  
  testthat::expect_equal(as.numeric(retval$edist1), 
                         2.881704,
                         tolerance = 1e-2)
  
  if (FALSE) {
  
  testthat::expect_equal(retval$soln2,
                         c(SBS1  = 1001.222,
                           SBS22 = 1999.372),
                         tolerance = 1) # Not sure why this needs to be so large
  
  testthat::expect_equal(as.numeric(retval$edist2),
                         1.414214,
                         tolerance = 1)
  }
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
                         3.176093,
                         tolerance = 1e-2)
  
  if (FALSE) {
  testthat::expect_equal(retval$soln2,
                         c(SBS5  = 1010.586, SBS10a = 2000.802),
                         tolerance = 1e-2)
  
  testthat::expect_equal(as.numeric(retval$edist2), 
                         2.897068,
                         tolerance = 1e-2)
  }
})


test_that("SparseAssignActivity1 (one spectrum) Test 3", {
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
  
  if (FALSE) {
  testthat::expect_equal(retval$soln2,
                         c(SBS3  = 297.3204,
                           SBS5  = 304.5565,
                           SBS4  = 301.2976,
                           SBS29 = 300.9959,
                           SBS24 = 296.2052,
                           SBS8  = 300.4491),
                         tolerance = 1)
  
  testthat::expect_equal(as.numeric(retval$edist2), 
                         1.732051,
                         tolerance = 1)
  }
})


test_that("SparseAssignActivity1 (one spectrum) Test 4", {
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
  
  if (FALSE) {
  testthat::expect_equal(retval$soln2,
                         c(SBS3  = 100.52494,
                           SBS5  = 100.22724,
                           SBS4  = 99.53495,
                           SBS29 = 94.38718,
                           SBS24 = 102.97684,
                           SBS8  = 102.51610),
                         tolerance = 1)
  
  testthat::expect_equal(as.numeric(retval$edist2), 2.854129,
                         tolerance = 1)
  }
})


test_that("SparseAssignActivity1 (one spectrum) Test 5", {
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
  
  if (FALSE) {
  testthat::expect_equal(retval$soln2,
                         c(SBS5  = 51.27930,
                           SBS29 = 66.66493,
                           SBS8  = 55.11371),
                         tolerance = 0.1)
  
  testthat::expect_equal(as.numeric(retval$edist2), 5.830952,
                         tolerance = 0.1)
  }
})

