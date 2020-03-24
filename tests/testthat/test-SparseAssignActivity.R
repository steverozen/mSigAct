context("test-SparseAssignActivity.R")

# Framework for testing \code{\link{SparseAssignActivity1}}.
# 
#' Tests only for single spectrum.
# 

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

