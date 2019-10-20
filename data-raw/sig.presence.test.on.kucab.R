# Test methyl eugenol


me.sig <- kucab.sigs[ , "Methyleugenol_1.25 mM", drop = FALSE]
# bad.control.sig <- kucab.sigs[ , "control", drop = FALSE]

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
  kucab.spectra[  , grep("Methyleugenol",
                         colnames(kucab.spectra),
                         fixed = TRUE)]

attr(me.spectra, "ref.genome") <- NULL
attr(me.sig, "ref.genome") <- NULL

m.opts <- DefaultManyOpts()

me1 <- me.spectra[ , 5, drop = FALSE] + 1e-100
fake.sig <- ICAMS::TransformCatalog(me1, target.ref.genome = "hg19",
                                    target.catalog.type = "counts.signature")

t1 <- 
  mSigAct::SignaturePresenceTest1(
    spectrum         = me1,
    sigs             = as.matrix(cbind(me.sig, kucab.control.sig)),
    # sigs              = as.matrix(cbind(me.sig, fake.sig)),
    target.sig.index = 1,
    m.opts           = m.opts,
    eval_f           = mSigAct:::ObjFnBinomMaxLH)
  