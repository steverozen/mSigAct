context("ReassignmentQP")

test_that("ReassignmentQP for SBS tumors", {
  spectra <- PCAWG7::spectra$PCAWG$SBS96
  exposure <- PCAWG7::exposure$PCAWG$SBS96
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  exp.QP <- 
    ReassignmentQP(spectra = spectra, 
                   exposure = exposure, 
                   sigs = sigs)
  
  expect_equal(sum(exp.QP[, 1] > 0), 5)
  expect_equal(sum(exposure[, 1] > 0), 5)
})