context("SigPresenceAssignActivity for SBS22")

test_that("SigPresenceAssignActivity for SBS Liver tumor Liver-HCC::SP97685", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  library(ICAMS)
  sbs.spectra <- PCAWG7::spectra$PCAWG$SBS96
  spectra <- sbs.spectra[, "Liver-HCC::SP97685", drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  
  sig.names <- setdiff(names(sigs.prop), cosmicsig::possible_artifacts())
  sig.names <- setdiff(sig.names, MSISigs()$SBS)
  
  
  output.dir <- "./tests/testthat/output/Liver-HCC.SP97685"
  dir.create(path = output.dir, recursive = TRUE)
  ICAMS::PlotCatalogToPdf(catalog = spectra,
                          file = file.path(output.dir, "Liver-HCC.SP97685.pdf"))
  
  sigs.to.use <- sigs[, sig.names, drop = FALSE]
  
  
  output.dir2 <- file.path(output.dir, "without.sbsh8")
  dir.create(path = output.dir2, recursive = TRUE)
  
  ICAMS::PlotCatalogToPdf(catalog = sigs.to.use,
                          file = file.path(output.dir2, "sbs.liver.signatures.pdf"))
  
  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 5
  my.opts$trace <- 1e15
  spaa.sbs22 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir2, "spaa.with.sbs22"),
                              max.level = ncol(sigs.to.use) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.sbs22 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir2, "sparse.with.sbs22"),
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
  
  pcawg.exp <- PCAWG7::exposure$PCAWG$SBS96
  pcawg.exp2 <- pcawg.exp[, "Liver-HCC::SP97685", drop = FALSE]
  pcawg.assign <- GetSampleSigActivity(spectra = spectra,
                                  exposure = pcawg.exp2,
                                  sigs = sigs, 
                                  sample.names = "Liver-HCC::SP97685", 
                                  output.dir = file.path(output.dir2, "pcawg.assign"))
  
  sigs.to.use2 <- sigs.to.use[, colnames(sigs.to.use) != "SBS22", drop = FALSE]
  spaa.no.sbs22 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use2,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir2, "spaa.without.sbs22"),
                              max.level = ncol(sigs.to.use2) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use2),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.no.sbs22 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use2,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir2, "sparse.no.sbs22"),
                         max.level = ncol(sigs.to.use2) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use2),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
  
  output.dir3 <- file.path(output.dir, "with.sbsh8")
  dir.create(path = output.dir3, recursive = TRUE)
  
  sbs.h8.file <-
    "~/packages/attribution_benchmarking/real_data/source_file/SBS_H8_SBS96.csv"
  sbs.catalog <- 
    ICAMS::ReadCatalog(file = sbs.h8.file, catalog.type = "counts.signature")
  sigs.to.use <- cbind(sigs.to.use, sbs.catalog)
  
  
  ICAMS::PlotCatalogToPdf(catalog = sigs.to.use,
                          file = file.path(output.dir3, "sbs.liver.signatures.pdf"))
  
  spaa.sbs22.sbsh8 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir3, "spaa.with.sbs22"),
                              max.level = ncol(sigs.to.use) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.sbs22.sbsh8 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir3, "sparse.with.sbs22"),
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
  
  
  sigs.to.use3 <- sigs.to.use[, colnames(sigs.to.use) != "SBS22", drop = FALSE]
  spaa.no.sbs22.sbsh8 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use3,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir3, "spaa.without.sbs22"),
                              max.level = ncol(sigs.to.use3) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use3),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.no.sbs22.sbsh8 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use3,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir3, "sparse.without.sbs22"),
                         max.level = ncol(sigs.to.use3) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use3),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
  
})

test_that("SigPresenceAssignActivity for SBS Liver tumor Liver-HCC::SP107101", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  library(ICAMS)
  sbs.spectra <- PCAWG7::spectra$PCAWG$SBS96
  spectra <- sbs.spectra[, "Liver-HCC::SP107101", drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  
  sig.names <- setdiff(names(sigs.prop), cosmicsig::possible_artifacts())
  sig.names <- setdiff(sig.names, MSISigs()$SBS)
  
  
  output.dir <- "./tests/testthat/output/Liver-HCC::SP107101"
  dir.create(path = output.dir, recursive = TRUE)
  ICAMS::PlotCatalogToPdf(catalog = spectra,
                          file = file.path(output.dir, "Liver-HCC::SP107101.pdf"))
  
  sigs.to.use <- sigs[, sig.names, drop = FALSE]
  
  
  output.dir2 <- file.path(output.dir, "without.sbsh8")
  dir.create(path = output.dir2, recursive = TRUE)
  
  ICAMS::PlotCatalogToPdf(catalog = sigs.to.use,
                          file = file.path(output.dir2, "sbs.liver.signatures.pdf"))
  
  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 5
  my.opts$trace <- 1e15
  spaa.sbs22 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir2, "spaa.with.sbs22"),
                              max.level = ncol(sigs.to.use) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.sbs22 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir2, "sparse.with.sbs22"),
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
  
  pcawg.exp <- PCAWG7::exposure$PCAWG$SBS96
  pcawg.exp2 <- pcawg.exp[, "Liver-HCC::SP107101", drop = FALSE]
  pcawg.assign <- GetSampleSigActivity(spectra = spectra,
                                       exposure = pcawg.exp2,
                                       sigs = sigs, 
                                       sample.names = "Liver-HCC::SP107101", 
                                       output.dir = file.path(output.dir2, "pcawg.assign"))
  
  sigs.to.use2 <- sigs.to.use[, colnames(sigs.to.use) != "SBS22", drop = FALSE]
  spaa.no.sbs22 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use2,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir2, "spaa.without.sbs22"),
                              max.level = ncol(sigs.to.use2) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use2),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.no.sbs22 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use2,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir2, "sparse.without.sbs22"),
                         max.level = ncol(sigs.to.use2) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use2),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
  
  output.dir3 <- file.path(output.dir, "with.sbsh8")
  dir.create(path = output.dir3, recursive = TRUE)
  
  sbs.h8.file <-
    "~/packages/attribution_benchmarking/real_data/source_file/SBS_H8_SBS96.csv"
  sbs.catalog <- 
    ICAMS::ReadCatalog(file = sbs.h8.file, catalog.type = "counts.signature")
  sigs.to.use <- cbind(sigs.to.use, sbs.catalog)
  
  ICAMS::PlotCatalogToPdf(catalog = sigs.to.use,
                          file = file.path(output.dir3, "sbs.liver.signatures.pdf"))
  
  spaa.sbs22.sbsh8 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir3, "spaa.with.sbs22"),
                              max.level = ncol(sigs.to.use) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.sbs22.sbsh8 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir3, "sparse.with.sbs22"),
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
  
  
  sigs.to.use3 <- sigs.to.use[, colnames(sigs.to.use) != "SBS22", drop = FALSE]
  spaa.no.sbs22.sbsh8 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use3,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir3, "spaa.without.sbs22"),
                              max.level = ncol(sigs.to.use3) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use3),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  sparse.no.sbs22.sbsh8 <-
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use3,
                         m.opts = my.opts,
                         output.dir = file.path(output.dir3, "sparse.without.sbs22"),
                         max.level = ncol(sigs.to.use3) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use3),
                         num.parallel.samples = 1,
                         mc.cores.per.sample = 180,
                         seed = 2561, 
                         max.subsets = .Machine$double.xmax)
})