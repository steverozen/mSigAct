context("test-AnySubsetPresent.R")

TestEsoSpectra <- function(indices = NULL) {
  eso.index <- grep("Eso", colnames(PCAWG7::PCAWG.WGS.SBS.96), fixed = TRUE)
  spectra <- PCAWG7::PCAWG.WGS.SBS.96[ , eso.index, drop = FALSE]
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
  return(PCAWG7::signature$genome$SBS96[ , sigs])
}

TestSignaturePresenceTestDouble <- function(extra.sig, eso.indices) {
  
  stopifnot(length(eso.indices) == 1)
  eso.spectra <- TestEsoSpectra(eso.indices)
  
  m.opts <- DefaultManyOpts()
  
  sigs.plus <- TestEsoSigs(extra.sig)
  set.seed(101010, kind = "L'Ecuyer-CMRG") 
  retval1 <- mSigAct::SignaturePresenceTest(
    spectra          = eso.spectra, 
    sigs             = sigs.plus,
    target.sig.index = 1, 
    m.opts           = m.opts, 
    eval_f           = ObjFnBinomMaxLHNoRoundOK, 
    mc.cores         = 1)
  
  set.seed(101010, kind = "L'Ecuyer-CMRG")
  retval2 <- SignaturePresenceTest1(
    spectrum         = eso.spectra,
    sigs             = sigs.plus,
    target.sig.index = 1,
    m.opts           = m.opts,
    eval_f           = ObjFnBinomMaxLHNoRoundOK)
  
  stopifnot(all.equal(
    retval1$`Eso-AdenoCA::SP111062`$chisq.p, 
    retval2$chisq.p))
  
  return(list(test = retval1, test1 = retval2))
}

TestAny1 <- function(extra.sig, eso.index) {
  
  eso.spectra <- TestEsoSpectra(eso.index)
  
  m.opts <- DefaultManyOpts()
  
  sigs.plus <- TestEsoSigs(extra.sig) # The extra signatures are signature names, and will be the first columns of sigs.plus
  
  set.seed(101010, kind = "L'Ecuyer-CMRG")  
  out <- AnySigSubsetPresent(spect           = eso.spectra,
                             all.sigs        = sigs.plus,
                             Ha.sigs.indices = 1:length(extra.sig),
                             eval_f          = mSigAct::ObjFnBinomMaxLHNoRoundOK,
                             m.opts          = m.opts)
  
  return(out)
}

test_that("TestAny1 and TestSignaturePresenceTest on SBS17a in esophageal sample 1", {
  any.retval <- TestAny1("SBS17a", 1)
  expect_equal(any.retval$all.Ha.info[[1]]$p, 0.09520268)
  # TestSigPresenceTests compares ressults for 
  spt.retval <- TestSignaturePresenceTestDouble("SBS17a", 1)
  expect_equal(spt.retval$test1$chisq.p, 0.09520268)
 
})

test_that("TestAny1 on SBS17a and SBS17b in esophageal sample 1", {
  any.retval2 <- TestAny1(c("SBS17a", "SBS17b"), 1)
  expect_equal(any.retval2$all.Ha.info[[3]]$sigs.added,
               "SBS17a,SBS17b")
  expect_equal(any.retval2$all.Ha.info[[1]]$sigs.added,
               "SBS17a")
  expect_equal(any.retval2$all.Ha.info[[1]]$p,
               0.09520268)
  expect_equal(any.retval2$all.Ha.info[[2]]$p,
               0.0009874463)
  expect_equal(any.retval2$all.Ha.info[[3]]$p,
               0.001310515)
})

# TODO STeve: test spectrum 6 and rounded reconstruction; is there a better way to decide not to use rounded reconstructon?