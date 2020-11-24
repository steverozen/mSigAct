
library(mSigAct)

mutation.type = "SBS96"

p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(
  PCAWG7::spectra$PCAWG[[mutation.type]])
lung <- p7$`Lung-SCC`

s96 <- PCAWG7::signature$genome$SBS96
l96 <- ncol(s96) - 1
s96 <- s96[ , -(l96:(l96+1))]

lung.subset <- lung[, c(1:4, 8, 18:22, 32:33)]

exposures <-
  apply(X          = lung.subset,
        MARGIN     = 2,
        FUN        = mSigAct::OptimizeExposureQP,
        signatures = s96)

ex01 <- exposures
ex01[ex01 >= 1] <- 1
ex01[ex01 < 1] <- 0
ex01.2 <- ex01[rowSums(ex01) > 0, ]
pdf(file = "black.white.pdf")
gplots::heatmap.2(ex01.2,
                  dendrogram = "none",
                  col = c("white", "black"),
                  Rowv = F, Colv = F, cexRow = 0.7,
                  trace = "none", key = F, labCol = F,
                  rowsep = 1:nrow(ex01.2))
dev.off()

signames1 <- rownames(exposures)[(which(rowSums(exposures) > 1000))]

ex1 <-
  apply(X          = lung.subset,
        MARGIN     = 2,
        FUN        = mSigAct::OptimizeExposureQP,
        signatures = s96[ , signames1])

ex2 <- ICAMSxtra::SortExposure(ex1)

ICAMS::PlotCatalogToPdf(
  ICAMS::as.catalog(
    lung.subset[, colnames(ex2)]), file = "lung.subset.pdf")

ICAMSxtra::PlotExposureToPdf(ex2, file = "bad.ex.labels.pdf")

colnames(ex2) <- NULL

ICAMSxtra::PlotExposureToPdf(ex2, file = "bad.ex.pdf")

