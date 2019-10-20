# Test methyl eugenol


me.sig <- kucab.sigs[ , "Methyleugenol_1.25 mM", drop = FALSE]
bad.control.sig <- kucab.sigs[ , "control", drop = FALSE]

kucab.control.idx <- grep("Control", colnames(kucab.spectra))
kucab.control.spectra <- kucab.spectra[ , kucab.control.idx]
kucab.control.sig <- rowMeans(kucab.control.spectra)
kucab.control.sig <- ICAMS::as.catalog(kucab.control.sig, 
                                       region = "genome",
                                       ref.genome = "hg19",
                                       catalog.type = "counts")
kucab.control.sig <- 
  ICAMS::TransformCatalog(kucab.control.sig,
                          target.catalog.type = "counts.signature")                                       

me.spectra <- 
  kucab.spectra[  , grep("Methyleugenol_(1.25 mM)",
                         colnames(kucab.spectra),
                         fixed = TRUE)]

attr(me.spectra, "ref.genome") <- NULL
attr(me.sig, "ref.genome") <- NULL

t1 <- 
  mSigAct::SignaturePresenceTest1(
    spectrum         = as.matrix(me.spectra[ , 1, drop = FALSE]),
    sigs             = as.matrix(cbind(me.sig, kucab.control.sig)),
    target.sig.index = 1,
    m.opts           = DefaultManyOpts(),
    eval_f           = mSigAct:::ObjFnBinomMaxLH)
