
BG.MEG.Test <- list()

BG.MEG.Test$sig.kucab <- mSigAct::kucab.sigs[ , "Methyleugenol..1.25.mM.", drop = FALSE]

# We will need this later, since the zou code eventually reports the
# extracted signature using the conventional row ordering used by ICAMS.
BG.MEG.Test$sig.ICAMS <- mSigAct::kucab.sigs.ICAMS[ , "Methyleugenol_1.25 mM", drop = FALSE]


me.spectra.kucab <- 
  mSigAct::kucab.sub.catalog[  , grep("MSM0.124", # Methyeugenol 1.25 mM
                                      colnames(mSigAct::kucab.sub.catalog),
                                      fixed = TRUE)]
rownames(me.spectra.kucab) <- mSigAct::kucab.sub.catalog[ ,1]

BG.MEG.Test$spectra.kucab <- me.spectra.kucab
rm(me.spectra.kucab)

BG.MEG.Test$x0.25.noise <- 
  mSigAct::Generate1KucabSynData(target.sig = BG.MEG.Test$sig.kucab,
                                 target.ratio.to.control = 0.25,
                                 num.replicates = 3, add.noise = TRUE)

BG.MEG.Test$x0.5.noise <- 
  mSigAct::Generate1KucabSynData(target.sig = BG.MEG.Test$sig.kucab,
                                 target.ratio.to.control = 0.5,
                                 num.replicates = 3, add.noise = TRUE)

BG.MEG.Test$x0.5.no.noise <- 
  mSigAct::Generate1KucabSynData(target.sig = BG.MEG.Test$sig.kucab,
                                 target.ratio.to.control = 0.5,
                                 num.replicates = 3, add.noise = FALSE)

usethis::use_data(BG.MEG.Test)
