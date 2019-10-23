context("SignaturePresenceTest for SBS96")

test_that("SignaturePresence1 96 1", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 0, tolerance = 1e-5)
})

test_that("SignaturePresence1 96 2", {
  input <- c(SBS5 = 10, SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 1, tolerance = 1e-5)
})

test_that("SignaturePresence1 96 3", {
  input <- c(SBS5 = 500, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 0.06987155, tolerance = 1e-5)
})

test_that("SignaturePresence1 96 4", {
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 0.007992491, tolerance = 1e-5)
})
