context("SignaturePresenceTest for SBS192")

all.abundance <- ICAMS::all.abundance

test_that("SignaturePresence1 SBS192 1", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0, tolerance = 1e-5)
})

test_that("SignaturePresence1 SBS192 2", {
  input <- c(SBS5 = 10, SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs = PCAWG7::signature$genome$SBS192)
  cat("test 2", retval$chisq.p, "\n")
  testthat::expect_equal(retval$chisq.p, 0.738, tolerance = 1e-3)
})

test_that("SignaturePresence1 SBS192 3", {
  input <- c(SBS5 = 500, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.5023189, tolerance = 1e-5)
})

test_that("SignaturePresence1 SBS192 4", {
  input <- c(SBS5 = 2000, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.03175961, tolerance = 1e-5)
})

test_that("SignaturePresence1 6", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.266199, tolerance = 1e-5)
})

test_that("SignaturePresence1 7", {
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000, SBS40 = 1000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs =  PCAWG7::signature$genome$SBS192)
  cat("test 7", retval$chisq.p, "\n")
  testthat::expect_equal(retval$chisq.p, 0.411, tolerance = 1e-2)
})

