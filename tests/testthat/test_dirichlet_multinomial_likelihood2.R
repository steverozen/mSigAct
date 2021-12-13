context("Use Dirichlet-multinomial distribution to calculate likelihood")

test_that("Use Dirichlet-multinomial distribution to calculate likelihood", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  prior.prop <- ExposureProportions("SBS96", cancer.type = "Lung-AdenoCA")
  prior.prop1 <- prior.prop[!names(prior.prop) %in% PossibleArtifacts()]
  sigs <- 
    ICAMS::ReadCatalog(file = "/home/e0012078/data/mSigAct_syn_data/ground.truth.syn.sigs.SBS96.csv",
                       catalog.type = "counts.signature")
  
  spectra <- ICAMS::ReadCatalog(file = "/home/e0012078/data/mSigAct_syn_data/ground.truth.syn.catalog.noisy.neg.binom.size.30.SBS96.csv")
  exposure <- ReadExposure(file = "/home/e0012078/data/mSigAct_syn_data/ground.truth.syn.exposures.noisy.neg.binom.size.30.SBS96.csv")
    
  sigs.to.use <- sigs[, names(prior.prop1), drop = FALSE]
  lung.tumor.indices <- grep(pattern = "Lung-AdenoCA", x = colnames(spectra))
  lung.spectra <- spectra[, lung.tumor.indices[2], drop = FALSE]
  lung.exposure <- RemoveZeroActivitySig(exposure[, lung.tumor.indices[2], drop = FALSE])
  
  # Use negative binomial distribution to calculate log likelihood
  sparse.out1 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = DefaultManyOpts(likelihood.dist = "neg.binom"),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), "lung.sparse.neg.binom"))
  
  expect_equal(sparse.out1$reconstruction.distances$sparse.assign.distances$cosine,
               0.9860534, tolerance = 1e-2)
  
  # Use multinomial distribution to calculate log likelihood (the default)
  sparse.out2 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), "lung.sparse.multinom"))
  
  expect_equal(sparse.out2$reconstruction.distances$sparse.assign.distances$cosine,
               0.9939246, tolerance = 1e-5)
  
  # Use Dirichlet-multinomial distribution to calculate log likelihood
  my.opts <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                             cp.factor = 100)
  sparse.out3 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 360,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.100"))
  
  expect_equal(sparse.out3$reconstruction.distances$sparse.assign.distances$cosine,
               0.9817432, tolerance = 1e-2)
  
  my.opts <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                             cp.factor = 1000)
  sparse.out4 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.1000"))
  
  expect_equal(sparse.out4$reconstruction.distances$sparse.assign.distances$cosine,
               0.9860423, tolerance = 1e-2)
  
  my.opts <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                             cp.factor = 10000)
  sparse.out5 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 306,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.10000"))
  
  expect_equal(sparse.out5$reconstruction.distances$sparse.assign.distances$cosine,
               0.9936892, tolerance = 1e-2)
  
  my.opts6 <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                             cp.factor = 100000)
  sparse.out6 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts6,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.100000"))
  
  my.opts7 <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                              cp.factor = 50000)
  sparse.out7 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts7,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.50000"))
  
  my.opts8 <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                              cp.factor = 75000)
  sparse.out8 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts8,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.75000"))
  
  my.opts9 <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                              cp.factor = 60000)
  sparse.out9 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts9,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.60000"))
  
  my.opts10 <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                              cp.factor = 65000)
  sparse.out10 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts10,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.65000"))
  
  my.opts11 <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                               cp.factor = 70000)
  sparse.out11 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts11,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.70000"))
  
  my.opts12 <- DefaultManyOpts(likelihood.dist = "dirichlet.multinom",
                               cp.factor = 75000)
  sparse.out12 <- 
    SparseAssignActivity(spectra = lung.spectra,
                         sigs = sigs.to.use,
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         m.opts = my.opts12,
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 60,
                         max.subsets = 1e15,
                         seed = 8386,
                         output.dir = file.path(tempdir(), 
                                                "lung.sparse.dirichlet.multinom.factor.75000"))
  
  expect_equal(sparse.out6$reconstruction.distances$sparse.assign.distances$cosine,
               c(0.9974084, 0.9935337), tolerance = 1e-2)
  
  unlink(file.path(tempdir(), "lung.sparse.neg.binom"), recursive = TRUE)
  unlink(file.path(tempdir(), "lung.sparse.multinom"), recursive = TRUE)
})