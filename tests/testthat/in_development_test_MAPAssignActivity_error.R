test_that("MAPAssignActivity for ID Catalog with errors -- not finished", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  catalog <- ICAMS::ReadCatalog(file = "testdata/PCAWG7-Prost-AdenoCA-ten-samples.csv")
  sample.index <- 1:2
  catID <- catalog[, sample.index, drop = FALSE]
  catID[ , 1] <- 0 # Try setting counts to 0 for 1 sample
  ID.sigs <- ICAMS::ReadCatalog(file = "testdata/COSMIC-v3-genome-ID-sigs.csv",
                                catalog.type = "counts.signature")
  mutation.type <- "ID"
  cancer.type <- "Prost-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- ID.sigs[, names(sigs.prop), drop = FALSE]
  output.dir <- file.path(tempdir(), "ID")
  
  retval <- MAPAssignActivity(
    spectra                 = catID,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = DefaultManyOpts(likelihood.dist = "neg.binom"),
    num.parallel.samples    = 1,
    mc.cores.per.sample     = 1,
    use.sparse.assign       = TRUE,
    use.sig.presence.test   = TRUE,
    seed                    = 8787)
  
  expect_true(TRUE)
 
  unlink(output.dir, recursive = TRUE)
  
})

test_that("Use sparse assignment for MAPAssignActivity", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  spectra <- PCAWG7::spectra$PCAWG$SBS96
  catalog.list <- PCAWG7::SplitPCAWGMatrixByTumorType(spectra)
  biliary.catalogs <- catalog.list$`Biliary-AdenoCA`
  
  sample.index <- 1
  catalog <- biliary.catalogs[, sample.index, drop = FALSE]
  SBS.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  mutation.type <- "SBS96"
  cancer.type <- "Biliary-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- SBS.sigs[, names(sigs.prop), drop = FALSE]
  my.opts <- DefaultManyOpts()
  
  output.dir <- file.path(tempdir(), paste0("test", 1:3))
  retval1 <- MAPAssignActivity(
    spectra                 = catalog,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir[1],
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 1,
    mc.cores.per.sample     = 30,
    seed                    = 2351)
  expect_equal(nrow(retval1$proposed.assignment), 14)
  
  retval2 <- MAPAssignActivity(
    spectra                 = catalog,
    sigs                    = sigs,
    output.dir              = output.dir[2],
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 1,
    mc.cores.per.sample     = 30,
    seed                    = 2351,
    use.sparse.assign       = TRUE)
  
  retval3 <- SparseAssignActivity(spectra = catalog,
                                  sigs = sigs,
                                  output.dir = output.dir[3],
                                  max.level = ncol(sigs) - 1,
                                  p.thresh = 0.01, 
                                  m.opts = my.opts, 
                                  num.parallel.samples = 1, 
                                  mc.cores.per.sample = 30,
                                  seed = 2351)
  expect_equal(retval2$proposed.assignment, retval3$proposed.assignment)
  
  sapply(output.dir, FUN = unlink, recursive = TRUE)
})