catSBS96.1 <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
my.sig.SBS96 <-
  ICAMS.shiny::CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                                     sig.type = "SBS96", region = "genome")
library(mSigAct)

mm <- DefaultManyOpts()
mm$trace <- 100

sigs.prop <- PCAWG7::exposure.stats$PCAWG$SBS96$`Biliary-AdenoCA`
sig.names <- rownames(sigs.prop)
sigs.prop <- unlist(sigs.prop[ , 2])
names(sigs.prop) <- sig.names


map.foo2 <- mSigAct::MAPAssignActivity1(spect = catSBS96.1, sigs = my.sig.SBS96,
                         sigs.presence.prop = sigs.prop,
                            max.level = 16,
                            p.thresh = 0.01,
                            eval_f = ObjFnBinomMaxLHNoRoundOK,
                            m.opts = mm #)
                            , max.mc.cores = 100 )# mc.cores.per.sample = 100)

small.sigs <-my.sig.SBS96[ , 1:7]
small.prop <- sigs.prop[colnames(small.sigs)]
small.prop[small.prop == 1] <-  0.99

# todo???

small.x <- mSigAct::MAPAssignActivity1(spect = catSBS96.1, sigs = small.sigs,
                                        sigs.presence.prop = small.prop,
                                        max.level = 4,
                                        p.thresh = 0.9,
                                        eval_f = ObjFnBinomMaxLHNoRoundOK,
                                        m.opts = mm #)
                                        , max.mc.cores = 100 )# mc.cores.per.sample = 100)
