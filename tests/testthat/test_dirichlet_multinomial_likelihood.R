context("Use Dirichlet-multinomial distribution to calculate likelihood")

test_that("Use Dirichlet-multinomial distribution to calculate likelihood", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  prior.prop <- ExposureProportions("SBS96", cancer.type = "Skin-Melanoma")
  prior.prop1 <- prior.prop[!names(prior.prop) %in% PossibleArtifacts()]
  sigs = PCAWG7::signature$genome$SBS96
  
  skin.spectra <- 
    PCAWG7::spectra$other.genome$SBS96[, c("Skin-Melanoma::ML_30_T_01",
                                           "Skin-Melanoma::ML_31_T_01"),drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.to.use <- sigs[, names(prior.prop1), drop = FALSE]
  
  # Use negative binomial distribution to calculate log likelihood
  sparse.out1 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         m.opts = DefaultManyOpts(likelihood.dist = "neg.binom"),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), "skin.sparse.neg.binom"))
  
  expect_equal(sparse.out1$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.8949538, 0.9083022), tolerance = 1e-2)
  
  # Use multinomial distribution to calculate log likelihood (the default)
  sparse.out2 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), "skin.sparse.multinom"))
  
  expect_equal(sparse.out2$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.917009, 0.994706), tolerance = 1e-5)
  
  DM.fit <- MGLM::MGLMfit(data = t(skin.spectra), dist = "DM")
  show(DM.fit)
  
  # Use Dirichlet-multinomial distribution to calculate log likelihood
  # sigs.to.use <- sigs.to.use[, 1:6]
  sparse.out3 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         m.opts = DefaultManyOpts(likelihood.dist = "dirichlet.multinom"),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), "skin.sparse.dirichlet.multinom"))
  
  expect_equal(sparse.out3$reconstruction.distances$sparse.assign.distances$cosine,
               0.9083022, tolerance = 1e-2)
  
  unlink(file.path(tempdir(), "skin.sparse.neg.binom"), recursive = TRUE)
  unlink(file.path(tempdir(), "skin.sparse.multinom"), recursive = TRUE)
})