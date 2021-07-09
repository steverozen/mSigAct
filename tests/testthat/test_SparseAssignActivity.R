context("SparseAssignActivity")

test_that("SparseAssignActivity for SBS Catalog", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
  sigs <- PCAWG7::signature$genome$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Lung-AdenoCA")
  sigs.to.use <- sigs[, names(sigs.prop), drop = FALSE]
  sparse.out <- SparseAssignActivity(spectra = spectra,
                                     sigs = sigs.to.use,
                                     output.dir = file.path(tempdir(), "Lung-AdenoCA"),
                                     max.level = ncol(sigs.to.use) - 1,
                                     p.thresh = 0.05 / ncol(spectra),
                                     num.parallel.samples = 2,
                                     mc.cores.per.sample = 30,
                                     seed = 2561)
  expect_equal(sparse.out$reconstruction.distances$sparse.assign.distances$cosine[1],
               0.9863608, tolerance = 1e-5)
  
  unlink(file.path(tempdir(), "Lung-AdenoCA"), recursive = TRUE)
  
})