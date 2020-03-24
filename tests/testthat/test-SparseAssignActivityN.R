context("test-SparseAssignActivityN.R")

SparseAssignTest <- function(sig.counts, trace = 0, 
                             mc.cores.per.sample = 1, 
                             num.parallel.samples = 1,
                             max.level = 5) {
  
  set.seed(101010, kind = "L'Ecuyer-CMRG")  
  
  if (!is.matrix(sig.counts)) {
    tmp.names <- names(sig.counts)
    sig.counts <- matrix(sig.counts, ncol = 1)
    rownames(sig.counts) <- tmp.names
  }
  sig.names <- rownames(sig.counts)
  colnames(sig.counts) <- paste0("tumor", 1:ncol(sig.counts))
  
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
  
  m.opts <- DefaultManyOpts()
  m.opts$trace <- trace
  
  SA.out <- 
    SparseAssignActivity(spectra              = spect,
                         sigs                 = some.sigs,
                         eval_f               = ObjFnBinomMaxLHMustRound,
                         m.opts               = m.opts,
                         mc.cores.per.sample  = mc.cores.per.sample,
                         num.parallel.samples = num.parallel.samples,
                         p.thresh = 0.5 # for backward testing compat.
    ) 
  return(SA.out)
}


test_that("SparseAssignActivity (multiple spectra) test 1", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <- SparseAssignTest(sig.counts = input, 
                             num.parallel.samples = 1,
                             mc.cores.per.sample  = 1) 
  expected <- matrix(c(999.2467, 1998.7533), ncol = 1)
  rownames(expected) <- names(input)
  colnames(expected) <- "tumor1"
  testthat::expect_equal(retval$exposure, expected, tolerance = 1e-2)
})


test_that("SparseAssignActivity (multiple spectra) test 2", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input.exp <- matrix(c(1000, 2000, 0, 2000, 10, 1000), ncol = 3)
  rownames(input.exp) <- c("SBS1", "SBS22")
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <-  SparseAssignTest(sig.counts = input.exp, 
                              num.parallel.samples = 1,
                              mc.cores.per.sample  = 1) 
  expected <- matrix(c(999.2467, 1998.7533, 0, 1996, 10.40123, 1000.59877),
                     ncol = 3)
  colnames(expected) <- paste0("tumor", 1:3)
  rownames(expected) <- rownames(input.exp)
  testthat::expect_equal(retval$exposure, expected, tolerance = 1e-2)
})

test_that("SparseAssignActivity (multiple spectra) test 3", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input <- c(SBS1 = 1000, SBS22 = 2000, SBS3 = 10, SBS40 = 10)
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <- SparseAssignTest(sig.counts = input, 
                             num.parallel.samples = 1,
                             mc.cores.per.sample  = 1) 
  expected <- matrix(c(999.2467, 1998.7533), ncol = 1)
  rownames(expected) <- names(input)
  colnames(expected) <- "tumor1"
  testthat::expect_equal(retval$exposure, expected, tolerance = 1e-2)
})

test_that("SparseAssignActivity (multiple spectra) test 4", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input <- c(SBS1 = 1000, SBS22 = 2000, SBS3 = 10, SBS40 = 10)
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <- SparseAssignTest(sig.counts = input, 
                             num.parallel.samples = 1,
                             mc.cores.per.sample  = 1,
                             max.level            = 1) 
  expected <- matrix(c(999.2467, 1998.7533), ncol = 1)
  rownames(expected) <- names(input)
  colnames(expected) <- "tumor1"
  testthat::expect_equal(retval$exposure, expected, tolerance = 1e-2)
})
