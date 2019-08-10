context("SparseAssignActivity")


test_that("SparseAssignTest1", {
  SparseAssignTest1 <- function() {
    retval <-  SparseAssignTestGeneric(
      sig.counts = c(SBS1 = 1000, SBS22 = 2000))
    testthat::expect_equal(retval$soln1,
                           c(SBS1 = 973.21485997560296, 
                             SBS22 = 2024.78514002439670),
                           tolerance = 10e-6)
    
    testthat::expect_equal(as.numeric(retval$edist1), 
                           15.25658687596772,
                           tolerance = 10e-6)
    
    testthat::expect_equal(retval$soln2,
                           c(SBS1  = 1001.0872211972701,
                             SBS22 = 1997.8671142054766),
                           tolerance = 1) # Not sure why this needs to be so large
    
    testthat::expect_equal(as.numeric(retval$edist2),
                           2.8248086460787443,
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
                           c(SBS3  = 0,
                             SBS5 = 1023.7736770381522, 
                             SBS10a = 1990.2263229618482),
                           tolerance = 10e-2)
    
    testthat::expect_equal(as.numeric(retval$edist1),
                           7.3353217488852822,
                           tolerance = 10e-2)
    
    testthat::expect_equal(retval$soln2,
                           c(SBS5  = 1017.8476499534260,
                             SBS10a = 2001.1238458342166),
                           tolerance = 10e-2)
    
    testthat::expect_equal(as.numeric(retval$edist2), 
                           3.0816408978467047,
                           tolerance = 10e-2)
    
    # return(retval)
  }
  SparseAssignTest2()
})


test_that("SparseAssignTest3", {
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
                             SBS8  = 409.47286575833471))
    
    testthat::expect_equal(as.numeric(retval$edist1), 29.058922589465876)
    
    testthat::expect_equal(retval$soln2,
                           c(SBS3  = 457.00823527475012,
                             SBS5  = 228.82984325255049,
                             SBS29 = 427.19460050811642,
                             SBS24 = 278.96098648348919,
                             SBS8  = 424.55158947103723))
    
    testthat::expect_equal(as.numeric(retval$edist2), 26.093099030576212)
    
    # return(retval)
  }
  SparseAssignTest3()
})


test_that("SparseAssignTest4", {
  skip("temp skip")
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
                             SBS8  = 117.97195281031038))
    
    testthat::expect_equal(as.numeric(retval$edist1), 9.6273257846362128)
    
    
    testthat::expect_equal(retval$soln2,
                           c(SBS5  = 140.36142697601778,
                             SBS4  = 147.98686390519100,
                             SBS24 = 162.52324641000195,
                             SBS8  = 139.74758200389428))
    
    testthat::expect_equal(as.numeric(retval$edist2), 9.1310602070161853)
    
    # return(retval)
  }
  SparseAssignTest4()
})

  
test_that("SparseAssignTest5", {
  skip("temp skip")
  SparseAssignTest5 <- function() {
    retval <- SparseAssignTestGeneric(
      sig.counts = c(SBS3=30, SBS5=30, SBS4=30, SBS29=30, SBS24=30, SBS8=30)
    )
    
    testthat::expect_equal(retval$soln1,
                           c(SBS3  = 0,
                             SBS5  = 56.488211683928093,
                             SBS4  = 0,
                             SBS29 = 70.317307606931379,
                             SBS24 = 0,
                             SBS8  = 50.194480709140528))
    
    testthat::expect_equal(as.numeric(retval$edist1), 5.6778792385700427)
    
    testthat::expect_equal(retval$soln2,
                           c(SBS5  = 51.289968752518106,
                             SBS29 = 66.599553081956088,
                             SBS8  = 55.241561633137891))
    
    testthat::expect_equal(as.numeric(retval$edist2), 5.6121626368613109)
    
    # return(retval)
  }
  SparseAssignTest5()
})

