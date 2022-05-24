context("GetSigActivity")

test_that("GetSigActivity for ID17 in PCAWG indel tumors", {
  spectra <- PCAWG7::spectra$PCAWG$ID
  exposure <- PCAWG7::exposure$PCAWG$ID
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
  sig.exp <- GetSigActivity(spectra = spectra, 
                            exposure = exposure,
                            sigs = sigs, 
                            sig.id = "ID17", 
                            output.dir = file.path(tempdir(), "ID17"))
  expect_equal(rowSums(sig.exp)["ID17"], 14, check.attributes = FALSE)
  unlink(file.path(tempdir(), "ID17"), recursive = TRUE)
})