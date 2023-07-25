test_that("MAPAssignActivity for ID Catalog", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  catalog <- ICAMS::ReadCatalog(file = "testdata/PCAWG7-Prost-AdenoCA-ten-samples.csv")
  sample.index <- 1
  catID <- catalog[, sample.index, drop = FALSE]
  ID.sigs <- ICAMS::ReadCatalog(file = "testdata/COSMIC-v3-genome-ID-sigs.csv",
                                catalog.type = "counts.signature")
  mutation.type <- "ID"
  cancer.type <- "Prost-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- ID.sigs[, names(sigs.prop), drop = FALSE]
  output.dir <- file.path(tempdir(), "ID")
  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 5
    
  retval <- MAPAssignActivity(
    spectra                 = catID,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 2,
    mc.cores.per.sample     = 30,
    seed                    = 8787)
  
  expect_true (all.equal(retval$proposed.assignment[, 1],
                         round(
                         c(ID1 = 60.0957928396838, 
                           ID2 = 12.5753318347774,  
                           ID3 = 20.5025441148954, 
                           ID5 = 74.6080111814531, 
                           ID8 = 19.7204896759284,  
                           ID9 = 15.131816682705, 
                           ID10 = 13.366013670557)),
                         tolerance = 1e-3))
  
  unlink(output.dir, recursive = TRUE)
})

test_that("Use sparse assignment for MAPAssignActivity", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  spectra <- PCAWG7::spectra$PCAWG$SBS96
  catalog.list <- PCAWG7::SplitPCAWGMatrixByTumorType(spectra)
  lung.catalogs <- catalog.list$`Lung-SCC`
  
  sample.index <- 1
  catalog <- lung.catalogs[, sample.index, drop = FALSE]
  SBS.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  mutation.type <- "SBS96"
  cancer.type <- "Lung-SCC"
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
    seed                    = 2351, 
    use.sparse.assign       = TRUE)
  expect_equal(nrow(retval1$proposed.assignment), 5)
  
  sapply(output.dir, FUN = unlink, recursive = TRUE)
})