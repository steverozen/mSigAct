context("Calculating distances between spectra and reconstructed spectra")

test_that("Identifying ID samples with low reconstruction accuracy", {
  ID.spectra <- PCAWG7::spectra$PCAWG$ID
  ID.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
  ID.exposure <- PCAWG7::exposure$PCAWG$ID
  
  indices <- 1:3
  distances <- CalculateDistance(spectra = ID.spectra[, indices],
                                 exposure = ID.exposure[, indices], 
                                 sigs = ID.sigs, 
                                 use.sparse.assign = TRUE)
  expect_equal(distances$scaled.manhattan[1], 0.3509820, tolerance = 1e-5)
})