catSBS96.1 <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
catSBS192.1 <- PCAWG7::spectra$PCAWG$SBS192[, 1, drop = FALSE]
my.exposure.SBS96 <- PCAWG7::exposure$PCAWG$SBS96[, 1, drop = FALSE]

catDBS78.1 <- PCAWG7::spectra$PCAWG$DBS78[, 1, drop = FALSE]
my.exposure.DBS78 <- PCAWG7::exposure$PCAWG$DBS78[, 1, drop = FALSE]

catID.1 <- PCAWG7::spectra$PCAWG$ID[, 1, drop = FALSE]
my.exposure.ID <- PCAWG7::exposure$PCAWG$ID[, 1, drop = FALSE]

my.sig.SBS96 <-
  CancerTypeToSigSubset(ca.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "SBS96", region = "genome")
my.sig.SBS192 <-
  CancerTypeToSigSubset(ca.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "SBS192", region = "genome")
my.sig.DBS78 <-
  CancerTypeToSigSubset(ca.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "DBS78", region = "genome")

my.sig.ID <-
  CancerTypeToSigSubset(ca.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "ID", region = "genome")


exposures.SBS96 <- GetExposureAndPlotToPdf(catalog = catSBS96.1,
                                     file = file.path(tempdir(), "test.SBS96.pdf"),
                                     sig.universe = my.sig.SBS96,
                                     num.of.replicates = 1000,
                                     method = decomposeQP)

sigfit.exposures.SBS96 <- GetExposureUseSigMiner(catalog = catSBS96.1,
                                           sig.universe = my.sig.SBS96)

exposures.SBS192 <- GetExposureAndPlotToPdf(catalog = catSBS192.1,
                                            file = file.path("x.pdf"), # tempdir(), "test.SBS192.pdf"),
                                            sig.universe = my.sig.SBS192,
                                            num.of.replicates = 1000,
                                            method = decomposeQP)
mysigs2 <- my.sig.SBS192[ , rownames(exposures.SBS192)]


foo3 <- mSigAct::SparseAssignActivity(spectra = catSBS192.1,
                                      sigs = mysigs2,
                                      max.level = 1,
                                      m.opts = mSigAct::DefaultManyOpts())

# foo3 <- mSigAct::AnySigSubsetPresent(spect = catSBS192.1, all.sigs = mysigs2,
#                                     Ha.sigs.indices = c(4, 6, 7, 8, 9, 10, 11),
#                                     m.opts = mSigAct::DefaultManyOpts())



sigfit.exposures.SBS192 <- GetExposureUseSigMiner(catalog = catSBS192.1,
                                                 sig.universe = my.sig.SBS192)

exposures.DBS78 <- GetExposureAndPlotToPdf(catalog = catDBS78.1,
                                           file = file.path(tempdir(), "test.DBS78.pdf"),
                                           sig.universe = my.sig.DBS78,
                                           num.of.replicates = 1000,
                                           method = decomposeQP)

sigfit.exposures.DBS78 <- GetExposureUseSigMiner(catalog = catDBS78.1,
                                                 sig.universe = my.sig.DBS78)

exposures.ID <- GetExposureAndPlotToPdf(catalog = catID.1,
                                           file = file.path(tempdir(), "test.ID.pdf"),
                                           sig.universe = my.sig.ID,
                                           num.of.replicates = 1000,
                                           method = decomposeQP)

sigfit.exposures.ID <- GetExposureUseSigMiner(catalog = catID.1,
                                                 sig.universe = my.sig.ID)

mo <- mSigAct::DefaultManyOpts()
mo$trace <- 1
mo$nbinom.size <- 3

library(profvis)
profvis(
foo5 <- mSigAct::SparseAssignActivity(spectra = catSBS192.1,
                                      sigs = mysigs2,
                                      max.level = 6,
                                      m.opts = mo))
