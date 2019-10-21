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


Generate1KucabSynDate <- function(target.sig, target.ratio.to.control, num.replicates, seed = NULL) {
  if (is.null(seed)) seed <- 10101
  controls <- sample(35, size = num.replicates, replace = TRUE)
  control.spectra <- kucab.control.spectra[ , controls, drop = FALSE]
  n.control.muts <- colSums(control.spectra)
  # START HERE
  
  
  
  
}
  