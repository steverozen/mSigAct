context("SparseAssignActivity")


test_that("SparseAssignTest1", {
  retval <-  SparseAssignTestGeneric(
    sig.counts = c(SBS1 = 1000, SBS22 = 2000), trace = 0)
  testthat::expect_equal(retval$soln1,
                         c(SBS1    = 999.2467, 
                           SBS22   = 1998.7533),
                         tolerance = 1e-2)
  
  testthat::expect_equal(as.numeric(retval$edist1), 
                         2.881704,
                         tolerance = 1e-2)
  
  testthat::expect_equal(retval$soln2,
                         c(SBS1  = 1001.222,
                           SBS22 = 1999.372),
                         tolerance = 1) # Not sure why this needs to be so large
  
  testthat::expect_equal(as.numeric(retval$edist2),
                         1.414214,
                         tolerance = 1)
})


test_that("SparseAssignTest2", {
  retval <- SparseAssignTestGeneric(
    sig.counts = c(SBS3 = 10, SBS5 = 1000, SBS10a = 2000),
    trace = 0)
  
  testthat::expect_equal(retval$soln1,
                         c(SBS3  = 0, SBS5 = 1011.794, SBS10a = 2002.206),
                         tolerance = 1e-2)
  
  testthat::expect_equal(as.numeric(retval$edist1),
                         3.176093,
                         tolerance = 1e-2)
  
  testthat::expect_equal(retval$soln2,
                         c(SBS5  = 1010.586, SBS10a = 2000.802),
                         tolerance = 1e-2)
  
  testthat::expect_equal(as.numeric(retval$edist2), 
                         2.897068,
                         tolerance = 1e-2)
})


test_that("SparseAssignTest3", {
  retval <- SparseAssignTestGeneric(
    sig.counts = 
      c(SBS3=300, SBS5=300, SBS4=300, SBS29=300, SBS24=300, SBS8=300),
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
})


test_that("SparseAssignTest4", {
  retval <- SparseAssignTestGeneric(
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
})


test_that("SparseAssignTest5", {
  retval <- SparseAssignTestGeneric(
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
  
  testthat::expect_equal(retval$soln2,
                         c(SBS5  = 51.27930,
                           SBS29 = 66.66493,
                           SBS8  = 55.11371),
                         tolerance = 0.1)
  
  testthat::expect_equal(as.numeric(retval$edist2), 5.830952,
                         tolerance = 0.1)
})


test_that("SparseAssignTest6", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <- XSparseAssignTestGeneric(sig.counts = input, mc.cores = 1)
  expected <- matrix(c(999.2467, 1998.7533), ncol = 1)
  rownames(expected) <- names(input)
  colnames(expected) <- "tumor1"
  testthat::expect_equal(retval, expected, tolerance = 1e-2)
})


test_that("SparseAssignTest7", {
  input.exp <- matrix(c(1000, 2000, 0, 2000, 10, 1000), ncol = 3)
  rownames(input.exp) <- c("SBS1", "SBS22")
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <-  XSparseAssignTestGeneric(sig.counts = input.exp, mc.cores = 1)
  expected <- matrix(c(999.2467, 1998.7533, 0, 1996, 10.40123, 1000.59877),
                     ncol = 3)
  colnames(expected) <- paste0("tumor", 1:3)
  rownames(expected) <- rownames(input.exp)
  testthat::expect_equal(retval, expected, tolerance = 1e-2)
})

