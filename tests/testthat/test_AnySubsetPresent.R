TestEsoSpectra <- function(indices = NULL) {
  ss <- PCAWG7::spectra$PCAWG$SBS96
  eso.index <- grep("Eso", colnames(ss), fixed = TRUE)
  spectra <- ss[ , eso.index, drop = FALSE]
  if (!is.null(indices)) {
    spectra <- spectra[ , indices, drop = FALSE]
  }
  return(spectra)
}

TestEsoSigs <- function(extra.sigs = NULL) {
  sigs <- c(
    "SBS1",
    "SBS2",
    "SBS3",
    "SBS5",
    "SBS13",
    "SBS18",
    "SBS28",
    "SBS40")
  if (!is.null(extra.sigs)) {
    sigs <- c(extra.sigs, sigs)
  }
  return(cosmicsig::COSMIC_v3.0$signature$GRCh37$SBS96[ , sigs])
}

TestSignaturePresenceTestDouble <- function(extra.sig, eso.indices) {

  stopifnot(length(eso.indices) == 1)
  eso.spectra <- TestEsoSpectra(eso.indices)

  m.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")

  sigs.plus <- TestEsoSigs(extra.sig)
  set.seed(101010, kind = "L'Ecuyer-CMRG")
  retval1 <- mSigAct::SignaturePresenceTest(
    spectra          = eso.spectra,
    sigs             = sigs.plus,
    target.sig.index = 1,
    m.opts           = m.opts,
    mc.cores         = 1)

  set.seed(101010, kind = "L'Ecuyer-CMRG")
  retval2 <- SignaturePresenceTest1(
    spectrum         = eso.spectra,
    sigs             = sigs.plus,
    target.sig.index = 1,
    m.opts           = m.opts)

  stopifnot(all.equal(
    retval1$`Eso-AdenoCA::SP111062`$chisq.p,
    retval2$chisq.p))

  return(list(test = retval1, test1 = retval2))
}

TestAny1 <- function(extra.sig, eso.index) {

  eso.spectra <- TestEsoSpectra(eso.index)

  m.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  m.opts$trace <- 100

  sigs.plus <- TestEsoSigs(extra.sig) # The extra signatures are signature names, and will be the first columns of sigs.plus

  set.seed(101010, kind = "L'Ecuyer-CMRG")
  out <- AnySigSubsetPresent(spect           = eso.spectra,
                             all.sigs        = sigs.plus,
                             Ha.sigs.indices = 1:length(extra.sig),
                             m.opts          = m.opts,
                             max.mc.cores    = 1) # Travis-CI will not use multiple cores

  return(out)
}

test_that("TestAny1 and TestSignaturePresenceTest on SBS17a in esophageal sample 1", {
  any.retval <- TestAny1("SBS17a", 1)
  expected <-0.0944294265780503
  expect_equal(any.retval$all.Ha.info[[1, "p"]], expected, tolerance = 1e-3)
  spt.retval <- TestSignaturePresenceTestDouble("SBS17a", 1)
  expect_equal(spt.retval$test1$chisq.p, expected, tolerance = 1e-3)
})


test_that("TestAny1 on SBS17a and SBS17b in esophageal sample 1", {
  any.retval2 <- TestAny1(c("SBS17a", "SBS17b"), 1)
  expect_equal(any.retval2$all.Ha.info[[3, "sigs.added"]],
               "SBS17a,SBS17b")
  expect_equal(any.retval2$all.Ha.info[[1, "sigs.added"]],
               "SBS17a")
  expect_equal(any.retval2$all.Ha.info[[1, "p"]],
               0.0944294265780503, tolerance = 1e-3)
  expect_equal(any.retval2$all.Ha.info[[2, "p"]],
               0.000910631136606777, tolerance = 1e-5)
  expect_equal(any.retval2$all.Ha.info[[3, "p"]],
               0.00139051152114872, tolerance = 1e-5)
})

