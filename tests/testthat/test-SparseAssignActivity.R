context("SparseAssignActivity")


test_that("SparseAssignTest1", {
  SparseAssignTest1 <- function() {
    retval <-  SparseAssignTestGeneric(
      sig.counts = c(SBS1 = 1000, SBS22 = 2000))
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
    
    # return(retval)
  }
  SparseAssignTest1()
})


test_that("SparseAssignTest2", {
  SparseAssignTest2 <- function() {
    retval <- SparseAssignTestGeneric(
      sig.counts = c(SBS3 = 10, SBS5 = 1000, SBS10a = 2000)
    )
    
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
                           tolerance = 10e-2)
    
    # return(retval)
  }
  SparseAssignTest2()
})


test_that("SparseAssignTest3", {
  skip("fix")
  SparseAssignTest3 <- function() {
    retval <- SparseAssignTestGeneric(
      sig.counts = c(SBS3=300, SBS5=300, SBS4=300, SBS29=300, SBS24=300, SBS8=300)
    )
    
    testthat::expect_equal(retval$soln1,
                           c(SBS3 = 340.97210938598397,
                             SBS5  = 278.90426791885159,
                             SBS4  = 0.00000000000000,
                             SBS29 = 418.20145869182011,
                             SBS24 = 350.44929824500957,
                             SBS8  = 409.47286575833471),
                           tolerance = 10e-2)
    
    testthat::expect_equal(as.numeric(retval$edist1),
                           29.058922589465876,
                           tolerance = 10e-2)
    
    testthat::expect_equal(retval$soln2,
                           c(SBS3  = 457.00823527475012,
                             SBS5  = 228.82984325255049,
                             SBS29 = 427.19460050811642,
                             SBS24 = 278.96098648348919,
                             SBS8  = 424.55158947103723),
                           tolerance = 1)
    
    testthat::expect_equal(as.numeric(retval$edist2), 
                           26.093099030576212,
                           tolerance = 1)
    
    # return(retval)
  }
  SparseAssignTest3()
})


test_that("SparseAssignTest4", {
  skip("fix")
  SparseAssignTest4 <- function() {
    retval <- SparseAssignTestGeneric(
      sig.counts = c(SBS3=100, SBS5=100, SBS4=100, SBS29=100, SBS24=100, SBS8=100)
    )
    
    testthat::expect_equal(retval$soln1,
                           c(SBS3  = 0,
                             SBS5  = 162.39808288294344,
                             SBS4  = 164.49837006361392,
                             SBS29 = 0,
                             SBS24 = 154.13159424313233,
                             SBS8  = 117.97195281031038),
                           tolerance = 10e-2)
    
    testthat::expect_equal(as.numeric(retval$edist1),
                           9.6273257846362128,
                           tolerance = 10e-2)
    
    
    testthat::expect_equal(retval$soln2,
                           c(SBS5  = 140.36142697601778,
                             SBS4  = 147.98686390519100,
                             SBS24 = 162.52324641000195,
                             SBS8  = 139.74758200389428),
                           tolerance = 1)
    
    testthat::expect_equal(as.numeric(retval$edist2),
                           9.1310602070161853,
                           tolerance = 1)
    
    # return(retval)
  }
  SparseAssignTest4()
})

  
test_that("SparseAssignTest5", {
  skip("fix")
  SparseAssignTest5 <- function() {
    retval <- SparseAssignTestGeneric(
      sig.counts = c(SBS3=30, SBS5=30, SBS4=30, SBS29=30, SBS24=30, SBS8=30),
      trace = 1
    )
    
    testthat::expect_equal(retval$soln1,
                           c(SBS3  = 0,
                             SBS5  = 56.488211683928093,
                             SBS4  = 0,
                             SBS29 = 70.317307606931379,
                             SBS24 = 0,
                             SBS8  = 50.194480709140528),
                           tolerance = 10e-2)
    
    testthat::expect_equal(as.numeric(retval$edist1),
                           5.6778792385700427,
                           tolerance = 10e-2)
    
    testthat::expect_equal(retval$soln2,
                           c(SBS5  = 51.289968752518106,
                             SBS29 = 66.599553081956088,
                             SBS8  = 55.241561633137891),
                           tolerance = 0.1)
    
    testthat::expect_equal(as.numeric(retval$edist2),
                           5.6121626368613109,
                           tolerance = 0.1)
    
    # return(retval)
  }
  SparseAssignTest5()
})


test_that("SparseAssignTest6", {
  skip("fix")
  input <- c(SBS1 = 1000, SBS22 = 2000)
  retval <-  XSparseAssignTestGeneric(sig.counts = input)
  expected <- matrix(c(973.21486, 2024.78514), ncol = 1)
  rownames(expected) <- names(input)
  colnames(expected) <- "tumor1"
  testthat::expect_equal(retval, expected, tolerance = 10e-6)
})


test_that("SparseAssignTest7", {
  skip("fix")
  input.exp <- matrix(c(1000, 2000, 0, 2000, 10, 1000), ncol = 3)
  rownames(input.exp) <- c("SBS1", "SBS22")
  retval <-  XSparseAssignTestGeneric(sig.counts = input.exp)
  expected <- matrix(c(973.2149, 2024.7851, 0, 1996, 8.745481, 1002.254519),
                     ncol = 3)
  colnames(expected) <- paste0("tumor", 1:3)
  rownames(expected) <- rownames(input.exp)
  testthat::expect_equal(retval, expected, tolerance = 10e-2)
})

