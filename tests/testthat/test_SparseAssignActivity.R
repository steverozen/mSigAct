context("SparseAssignActivity")

test_that("SparseAssignActivity for SBS Catalog", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
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
  
  MAP.out <- MAPAssignActivity(spectra = spectra,
                               sigs = sigs.to.use,
                               output.dir = file.path(tempdir(), "Lung-AdenoCA2"),
                               max.level = ncol(sigs.to.use) - 1,
                               p.thresh = 0.05 / ncol(spectra),
                               num.parallel.samples = 2,
                               mc.cores.per.sample = 30,
                               seed = 2561,
                               use.sparse.assign = TRUE)
  expect_equal(sparse.out$reconstruction.distances$sparse.assign.distances$cosine[1],
               0.9863608, tolerance = 1e-5)
  expect_equal(sparse.out$proposed.assignment, MAP.out$proposed.assignment)
  
  unlink(file.path(tempdir(), "Lung-AdenoCA"), recursive = TRUE)
  unlink(file.path(tempdir(), "Lung-AdenoCA2"), recursive = TRUE)
})

test_that("SparseAssignActivity for DBS Catalog", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$DBS78))
  spectra <- PCAWG7::spectra$PCAWG$DBS78[, indices[1], drop = FALSE]
  sigs <- PCAWG7::signature$genome$DBS78
  sigs.prop <- ExposureProportions(mutation.type = "DBS78",
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
  
  MAP.out <- MAPAssignActivity(spectra = spectra,
                               sigs = sigs.to.use,
                               output.dir = file.path(tempdir(), "Lung-AdenoCA2"),
                               max.level = ncol(sigs.to.use) - 1,
                               p.thresh = 0.05 / ncol(spectra),
                               num.parallel.samples = 2,
                               mc.cores.per.sample = 30,
                               seed = 2561,
                               use.sparse.assign = TRUE)
  expect_equal(sparse.out$reconstruction.distances$sparse.assign.distances$cosine[1],
               0.9992178, tolerance = 1e-5)
  expect_equal(sparse.out$proposed.assignment, MAP.out$proposed.assignment)
  
  unlink(file.path(tempdir(), "Lung-AdenoCA"), recursive = TRUE)
  unlink(file.path(tempdir(), "Lung-AdenoCA2"), recursive = TRUE)
})

test_that("SparseAssignActivity for ID Catalog", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$ID))
  spectra <- PCAWG7::spectra$PCAWG$ID[, indices[1], drop = FALSE]
  sigs <- PCAWG7::signature$genome$ID
  sigs.prop <- ExposureProportions(mutation.type = "ID",
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
  
  MAP.out <- MAPAssignActivity(spectra = spectra,
                               sigs = sigs.to.use,
                               output.dir = file.path(tempdir(), "Lung-AdenoCA2"),
                               max.level = ncol(sigs.to.use) - 1,
                               p.thresh = 0.05 / ncol(spectra),
                               num.parallel.samples = 2,
                               mc.cores.per.sample = 30,
                               seed = 2561,
                               use.sparse.assign = TRUE)
  expect_equal(sparse.out$reconstruction.distances$sparse.assign.distances$cosine[1],
               0.9780138, tolerance = 1e-5)
  expect_equal(sparse.out$proposed.assignment, MAP.out$proposed.assignment)
  
  unlink(file.path(tempdir(), "Lung-AdenoCA"), recursive = TRUE)
  unlink(file.path(tempdir(), "Lung-AdenoCA2"), recursive = TRUE)
})