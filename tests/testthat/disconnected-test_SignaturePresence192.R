m.opts <- DefaultManyOpts()

test_that("SignaturePresence1 SBS192 1", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input,
    input.sigs = PCAWG7::signature$genome$SBS192,
    m.opts     = m.opts)
  testthat::expect_equal(retval$chisq.p, 0, tolerance = 1e-5)
})

test_that("SignaturePresence1 SBS192 2", {
  input <- c(SBS5 = 10, SBS1 = 1000, SBS22 = 2000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input,
    input.sigs = PCAWG7::signature$genome$SBS192,
    m.opts     = m.opts)
  testthat::expect_equal(retval$chisq.p, 0.7404775, tolerance = 1e-5)
})

test_that("SignaturePresence1 SBS192 3", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input <- c(SBS5 = 500, SBS3 = 10000, SBS22 = 20000)
  set.seed(101010, kind = "L'Ecuyer-CMRG")
  retval <- TestSignaturePresenceTest1(
    sig.counts = input,
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p,
                         0.4904425,
                         tolerance = 1e-5)
})

test_that("SignaturePresence1 SBS192 4", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input <- c(SBS5 = 2000, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input,
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.03174639, tolerance = 1e-4)
})

test_that("SignaturePresence1 SBS192 5", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input,
    input.sigs = PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.266199, tolerance = 1e-3)
})

test_that("SignaturePresence1 SBS192 6", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input <- c(SBS5 = 800, SBS3 = 10000, SBS22 = 20000, SBS40 = 1000)
  retval <- TestSignaturePresenceTest1(
    sig.counts = input,
    input.sigs =  PCAWG7::signature$genome$SBS192)
  testthat::expect_equal(retval$chisq.p, 0.3901646, # 0.4110288,
                         tolerance = 1e-4)
  testthat::expect_equal(retval$with, -905.5161, tolerance = 0.01)
  testthat::expect_equal(retval$without, -905.8853, tolerance = 0.01)
})

