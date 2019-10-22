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

colnames(kucab.control.sig) <- "control.average.sig"

me.spectra <- 
  kucab.spectra[  , grep("Methyleugenol",
                         colnames(kucab.spectra),
                         fixed = TRUE)]

attr(me.spectra, "ref.genome") <- NULL
attr(me.sig, "ref.genome") <- NULL

m.opts <- DefaultManyOpts()

fake.sig <- ICAMS::TransformCatalog(me1, target.ref.genome = "hg19",
                                    target.catalog.type = "counts.signature")

sbs5 <- PCAWG7::signature$genome$SBS96[ , "SBS5"]

me.sig0 <- me.sig
colnames(me.sig0) <- "MEG.sig0"
me.sig0[17:32, ] <- 0
me.sig0[49:96, ] <- 0
ICAMS::PlotCatalog(me.sig0)

me.sig1 <- me.sig
colnames(me.sig1) <- "MEG.sig1"
me.sig1[49:64, ] <- 0
ICAMS::PlotCatalog(me.sig1)

flat.sig <- as.matrix(rep(1/96, 96))
colnames(flat.sig) <- "flat.sig"

me1 <- me.spectra[ , 3, drop = FALSE] + 1e-100
t1 <- 
  mSigAct::SignaturePresenceTest1(
    spectrum         = me1[1:16, ],
    sigs             = as.matrix(cbind(me.sig, kucab.control.sig))[1:16, ],
    # sigs              = as.matrix(cbind(me.sig, fake.sig)),
    target.sig.index = 1,
    m.opts           = m.opts,
    eval_f           = mSigAct:::ObjFnBinomMaxLH)


Generate1KucabSynDate <- 
  function(target.sig, target.ratio.to.control, num.replicates, seed = NULL) {
  if (is.null(seed)) seed <- 10101
  set.seed(seed)
  controls <- sample(35, size = num.replicates, replace = TRUE)
  control.spectra <- kucab.control.spectra[ , controls, drop = FALSE]

  num.control.muts <- colSums(control.spectra)
  
  num.target.muts <- num.control.muts * target.ratio.to.control
  
  target.spectra.no.noise <- target.sig %*% matrix(num.target.muts, nrow = 1)
  
  target.spectra.noise <- 
    matrix(sapply(
      round(target.spectra.no.noise),
      function(mu) rnbinom(n = 1, size = 5, mu = mu)),
      nrow = nrow(target.spectra.no.noise))
  
  test.spectra <- control.spectra + target.spectra.noise
  
  colnames(test.spectra) <- paste0("test.spec.", 1:ncol(test.spectra))
  
  # Need to return the exposures too
  return(list(test.spectra = test.spectra,
              bg.exposures = num.control.muts,
              target.exposures = colSums(target.spectra.noise)))
}
  
SPIKE.TEST <- Generate1KucabSynDate(target.sig = me.sig,
                             target.ratio.to.control = 0.1,
                             num.replicates = 3)  
  
  
tt1.1 <- 
  mSigAct::SignaturePresenceTest1(
    spectrum         = SPIKE.TEST$test.spectra[ ,1, drop = FALSE],
    sigs             = as.matrix(cbind(me.sig, kucab.control.sig)),
    target.sig.index = 1,
    m.opts           = m.opts,
    eval_f           = mSigAct:::ObjFnBinomMaxLH)  
  

foo.0.1 <- FindSignatureMinusBackground(SPIKE.TEST$test.spectra, 
           mSigAct::kucab.control.bg,
           algorithm = 'NLOPT_LN_COBYLA',
           # algorithm = "NLOPT_GN_DIRECT", # does not work well
           maxeval = 10000, 
           print_level = 0,
           xtol_rel = 0.001,  # 0.0001,)
           xtol_abs = 0.0001,
           start.b.fraction = 0.67)

foo.0.1.p <- mSigAct::SignaturePresenceTest1(
     spectrum = SPIKE.TEST$test.spectra[ , 1, drop = FALSE],
     sigs    = as.matrix(cbind(mSigAct::kucab.control.bg$background.sig, foo.0.1$target.sig)),
     target.sig.index = 1,
     m.opts = m.opts,
     eval_f = mSigAct:::ObjFnBinomMaxLH)
