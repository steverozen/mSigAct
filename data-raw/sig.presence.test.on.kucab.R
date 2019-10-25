# Test methyleugenol

if (FALSE) { # This signature is already stored in mSigAct::kucab.control.bg$background.sig

kucab.control.idx <- grep("Control", colnames(kucab.spectra.ICAMS))
kucab.control.spectra.ICAMS <- kucab.spectra.ICAMS[ , kucab.control.idx]

kucab.control.sig.ICAMS <- rowMeans(kucab.control.spectra.ICAMS)
kucab.control.sig.ICAMS <- ICAMS::as.catalog(kucab.control.sig.ICAMS, 
                                       region = "genome",
                                       catalog.type = "counts")
kucab.control.sig.ICAMS <- 
  ICAMS::TransformCatalog(kucab.control.sig.ICAMS,
                          target.catalog.type = "counts.signature")

colnames(kucab.control.sig.ICAMS) <- "control.average.sig"
}

me.sig.kucab <- mSigAct::kucab.sigs[ , "Methyleugenol..1.25.mM.", drop = FALSE]

me.spectra.kucab <- 
  mSigAct::kucab.sub.catalog[  , grep("MSM0.124", # Methyeugenol 1.25 mM
                         colnames(mSigAct::kucab.sub.catalog),
                         fixed = TRUE)]
rownames(me.spectra.kucab) <- mSigAct::kucab.sub.catalog[ ,1]

kucab.control.spectra.kucab <- 
  mSigAct::kucab.sub.catalog[ , as.character(mSigAct::kucab.muts.control$Sample)]



# me.sig0 <- me.sig
# colnames(me.sig0) <- "MEG.sig0"
# me.sig0[17:32, ] <- 0
# me.sig0[49:96, ] <- 0
# ICAMS::PlotCatalog(me.sig0)

# me.sig1 <- me.sig
# colnames(me.sig1) <- "MEG.sig1"
# me.sig1[49:64, ] <- 0
# ICAMS::PlotCatalog(me.sig1)

# flat.sig <- as.matrix(rep(1/96, 96))
# colnames(flat.sig) <- "flat.sig"

Generate1KucabSynData <- 
  function(target.sig, target.ratio.to.control, num.replicates, seed = NULL) {
  if (is.null(seed)) seed <- 10101
  set.seed(seed)
  controls <- sample(35, size = num.replicates, replace = TRUE)
  control.spectra <- kucab.control.spectra.kucab[ , controls, drop = FALSE]

  num.control.muts <- colSums(control.spectra)
  
  num.target.muts <- num.control.muts * target.ratio.to.control
  
  target.spectra.no.noise <- target.sig %*% matrix(num.target.muts, nrow = 1)
  
  target.spectra.noise <- 
    matrix(sapply(
      round(target.spectra.no.noise),
      function(mu) rnbinom(n = 1, size = 10000, mu = mu)),
      nrow = nrow(target.spectra.no.noise))
  
  test.spectra <- control.spectra + target.spectra.noise
  
  colnames(test.spectra) <- paste0("test.spec.", 1:ncol(test.spectra))
  
  # Need to return the exposures too
  return(list(test.spectra = test.spectra,
              bg.exposures = num.control.muts,
              target.exposures = colSums(target.spectra.noise)))
}
  
SPIKE.TEST <- Generate1KucabSynData(target.sig = me.sig,
                             target.ratio.to.control = 0.5,
                             num.replicates = 3)
  
  

# mSigAct tests below

KucabToICAMSSpectra <- function(m) {
  stopifnot(rownames(m) == rownames(kucab.sub.catalog))
  stopifnot(is.numeric(m[ ,1]))
  new.rownames <- ICAMS:::Unstaple96(rownames(m))
  rownames(m) <- new.rownames
  m <- m[ICAMS::catalog.row.order$SBS96, ]
  return(ICAMS::as.catalog(m, region = "genome", catalog.type = "counts"))
}

test.spectra <- KucabToICAMSSpectra(SPIKE.TEST$test.spectra)
foo.0.1 <- FindSignatureMinusBackground(
  test.spectra, 
  mSigAct::kucab.control.bg,
  algorithm = 'NLOPT_LN_COBYLA',
  # algorithm = "NLOPT_GN_DIRECT", # does not work well
  maxeval = 10000, 
  print_level = 0,
  xtol_rel = 0.001,  # 0.0001,)
  xtol_abs = 0.0001,
  start.b.fraction = 0.67)

me.sig.ICAMS <- kucab.sigs.ICAMS[ , "Methyleugenol_1.25 mM", drop = FALSE]


m.opts <- DefaultManyOpts()
foo.0.1.p <- mSigAct::SignaturePresenceTest1(
     spectrum = test.spectra[ , 1, drop = FALSE],
     sigs    = as.matrix(cbind(mSigAct::kucab.control.bg$background.sig, foo.0.1$target.sig)),
     target.sig.index = 1,
     m.opts = m.opts,
     eval_f = mSigAct:::ObjFnBinomMaxLH)



t1 <- 
  mSigAct::SignaturePresenceTest1(
    spectrum         = me1[1:16, ],
    sigs             = as.matrix(cbind(me.sig, kucab.control.sig))[1:16, ],
    # sigs              = as.matrix(cbind(me.sig, fake.sig)),
    target.sig.index = 1,
    m.opts           = m.opts,
    eval_f           = mSigAct:::ObjFnBinomMaxLH)



tt1.1 <- 
  mSigAct::SignaturePresenceTest1(
    spectrum         = SPIKE.TEST$test.spectra[ ,1, drop = FALSE],
    sigs             = as.matrix(cbind(me.sig, kucab.control.sig)),
    target.sig.index = 1,
    m.opts           = m.opts,
    eval_f           = mSigAct:::ObjFnBinomMaxLH)  

