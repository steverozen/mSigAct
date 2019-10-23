context("SignaturePresenceTest for SBS192")

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
  testthat::expect_equal(retval$chisq.p, 1, tolerance = 1e-5)
})

test_that("SignaturePresence1 SBS192 3", {
  input <- c(SBS5 = 500, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.5503635, tolerance = 1e-5)
})

test_that("SignaturePresence1 SBS192 4", {
  input <- c(SBS5 = 2000, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input, 
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.02880872, tolerance = 1e-5)
})

if (FALSE) {
test_that("SignaturePresence1 6", {
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(sig.counts = input, 
                                       input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.266199, tolerance = 1e-5)
})

test_that("SignaturePresence1 7", {
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000, SBS40 = 1000)
  retval <- TestSignaturePresenceTest1(sig.counts = input, 
                                       input.sigs =  PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.266199, tolerance = 1e-5)
})
}
