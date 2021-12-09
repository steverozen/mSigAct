context("Use Dirichlet-multinomial distribution to calculate likelihood")

test_that("Use Dirichlet-multinomial distribution to calculate likelihood", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  prior.prop <- ExposureProportions("SBS96", cancer.type = "Skin-Melanoma")
  prior.prop1 <- prior.prop[!names(prior.prop) %in% PossibleArtifacts()]
  sigs <- 
    ICAMS::ReadCatalog(file = "/home/e0012078/data/mSigAct_syn_data/ground.truth.syn.sigs.SBS96.csv",
                       catalog.type = "counts.signature")
  
  spectra <- ICAMS::ReadCatalog(file = "/home/e0012078/data/mSigAct_syn_data/ground.truth.syn.catalog.noisy.neg.binom.size.30.SBS96.csv")
  exposure <- ReadExposure(file = "/home/e0012078/data/mSigAct_syn_data/ground.truth.syn.exposures.noisy.neg.binom.size.30.SBS96.csv")
    
  sigs.to.use <- sigs[, names(prior.prop1), drop = FALSE]
  skin.tumor.indices <- grep(pattern = "Skin", x = colnames(spectra))
  skin.spectra <- spectra[, skin.tumor.indices[1:2]]
  skin.exposure <- RemoveZeroActivitySig(exposure[, skin.tumor.indices[1:2]])
  
  # Use negative binomial distribution to calculate log likelihood
  sparse.out1 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = DefaultManyOpts(likelihood.dist = "neg.binom"),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), "skin.sparse.neg.binom"))
  
  expect_equal(sparse.out1$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.9965910, 0.9860534), tolerance = 1e-2)
  
  # Use multinomial distribution to calculate log likelihood (the default)
  sparse.out2 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), "skin.sparse.multinom"))
  
  expect_equal(sparse.out2$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.9974528, 0.9939246), tolerance = 1e-5)
  
  # Use Dirichlet-multinomial distribution to calculate log likelihood
  sparse.out3 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = DefaultManyOpts(likelihood.dist = "dirichlet.multinom"),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "skin.sparse.dirichlet.multinom.factor.100"))
  
  expect_equal(sparse.out3$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.9661144, 0.9817432), tolerance = 1e-2)
  
  sparse.out4 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = DefaultManyOpts(likelihood.dist = "dirichlet.multinom"),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "skin.sparse.dirichlet.multinom.factor.200"))
  
  expect_equal(sparse.out4$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.9899887, 0.9860423), tolerance = 1e-2)
  
  sparse.out5 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = DefaultManyOpts(likelihood.dist = "dirichlet.multinom"),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "skin.sparse.dirichlet.multinom.factor.1000"))
  
  expect_equal(sparse.out5$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.9973036, 0.9936892), tolerance = 1e-2)
  
  sparse.out6 <- 
    SparseAssignActivity(spectra = skin.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = DefaultManyOpts(likelihood.dist = "dirichlet.multinom"),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "skin.sparse.dirichlet.multinom.factor.1500"))
  
  expect_equal(sparse.out6$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.9974084, 0.9935337), tolerance = 1e-2)
  
  unlink(file.path(tempdir(), "skin.sparse.neg.binom"), recursive = TRUE)
  unlink(file.path(tempdir(), "skin.sparse.multinom"), recursive = TRUE)
})