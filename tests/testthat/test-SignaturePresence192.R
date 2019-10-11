context("SignaturePresenceTest")

test_that("SignaturePresence1 1", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(sig.counts = input, 
                                       input.sigs=PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0, tolerance = 1e-5)
})

test_that("SignaturePresence1 2", {
  input <- c(SBS5 = 10, SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 1, tolerance = 1e-5)
})

test_that("SignaturePresence1 3", {
  input <- c(SBS5 = 13, SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 0.7847336, tolerance = 1e-5)
})

test_that("SignaturePresence1 4", {
  input <- c(SBS5 = 20, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 0.7617834, tolerance = 1e-5)
})

test_that("SignaturePresence1 5", {
  input <- c(SBS5 = 500, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 0.08645732, tolerance = 1e-5)
})

test_that("SignaturePresence1 6", {
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(sig.counts = input)
  testthat::expect_equal(retval$chisq.p, 0.009854126, tolerance = 1e-5)
})
