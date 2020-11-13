p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(PCAWG7::spectra$PCAWG$SBS96)
cancer.type <- "Liver-HCC"
lungscc <- p7[[cancer.type]]

library(mSigAct)

mm <- DefaultManyOpts()
mm$trace <- 100

sigs.prop <- PCAWG7::exposure.stats$PCAWG$SBS96[[cancer.type]]
sig.names <- rownames(sigs.prop)
sigs.prop <- unlist(sigs.prop[ , 2])
names(sigs.prop) <- sig.names

spect <- lungscc[ , 1, drop = FALSE]

mapout <-
  mSigAct::MAPAssignActivity1(
    spect = spect,
    sigs = PCAWG7::signature$genome$SBS96[ , names(sigs.prop)],
    sigs.presence.prop = sigs.prop,
    max.level = length(sigs.prop) - 1,
    p.thresh = 0.01,
    eval_f = ObjFnBinomMaxLHNoRoundOK,
    m.opts = mm #)
    , max.mc.cores = 100 )# mc.cores.per.sample = 100)
# todo remove p value and other non-useful columns from output


xx <- ListOfList2Tibble(mapout$everything)

dplyr::arrange(xx, MAP)[(nrow(xx) - 5):nrow(xx), ]

# select.best and todo compare to PCAWG exxxposure

best <- dplyr::arrange(xx, MAP)[nrow(xx),  ]

px <- PCAWG7::exposure$PCAWG$SBS96[ , colnames(spect)]
px <- px[px > 0]
pxt <- tibble(sig.id = names(px), px)
bx <- best$exp[[1]]
bxt <- tibble(sig.id = names(bx), bx )
most.sparse.exp <- xx[175, "exp"][[1]][[1]]
most.sparse <- tibble(sig.id = names(most.sparse.exp), most.sparse = most.sparse.exp)

QP.exp <- OptimizeExposureQP(spect, PCAWG7::signature$genome$SBS96[ , names(bx)])
qpt <- tibble(sig.id = names(QP.exp), qp = QP.exp)

QP.sparse.exp <- OptimizeExposureQP(spect, PCAWG7::signature$genome$SBS96[ , names(most.sparse.exp)])
qp.sparse.t <- tibble(sig.id = names(QP.sparse.exp), qp.sparse = QP.sparse.exp)

library(dplyr)
left_join(left_join(full_join(right_join(pxt, bxt), qpt), most.sparse), qp.sparse.t)



r.p <- ReconstructSpectrum(PCAWG7::signature$genome$SBS96, exp = px, use.sig.names = TRUE) # PCAWG
r.b <- ReconstructSpectrum(PCAWG7::signature$genome$SBS96, exp = bx, use.sig.names = TRUE) # MAP best
# Minimal number of signatures

r.qp <- ReconstructSpectrum(PCAWG7::signature$genome$SBS96, exp = qp[ ,1], use.sig.names = TRUE)
minb.b <- ReconstructSpectrum(PCAWG7::signature$genome$SBS96, exp = most.sparse.exp, use.sig.names = TRUE) # MAP most sparse
r.most.sparse <- ReconstructSpectrum(PCAWG7::signature$genome$SBS96, exp = QP.sparse.exp, use.sig.names = TRUE) # MAP most sparse


sol.matrix <- cbind(spect, r.p, r.b, r.qp, minb.b, r.most.sparse)
colnames(sol.matrix) <- c("spect", "PCAWG7", "MAP.best",  "MAP.best.QP.obj.fn", "sparse", "sparse.QP.obj.fn")

ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(sol.matrix), file = "foo.pdf")

philentropy::distance(t(sol.matrix), use.row.names = TRUE)

philentropy::distance(t(sol.matrix), use.row.names = TRUE, method = "cosine")

# to do reoptimize minimal solution using cossim or euclidean distance



