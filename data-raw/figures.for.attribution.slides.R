# This strictly "as-is" code for figures showing
#
# 1. attribution using all signatures
#
# 2. the above, but with signatures filtered for a minium number of
#    attributed mutations
#
# 3. attribution using only the signatures previously found in
#    a particular cancer type
#
# 4. attribution based on searching all subsets of signatures
#    that could plausibly explain the given spectrum, followed
#    by selection according to maximum a posteriori probabiliy

library(mSigAct)
library(PCAWG7)
library(ICAMS)
library(ICAMSxtra)
library(tibble)
library(dplyr)

mutation.type = "SBS96"

p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(
  PCAWG7::spectra$PCAWG[[mutation.type]])

# We will be using squamous cell lung cancers as examples
lung <- p7$`Lung-SCC`

s96 <- PCAWG7::signature$genome$SBS96
l96 <- ncol(s96) - 1
s96 <- s96[ , -(l96:(l96+1))]

# A subset of the Lung-SCC with similar mutation load for
# ease plotting
lung.subset <- lung[, c(1:4, 8, 18:22, 32:33)]

# 1. Attribute using all signatures
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

# 2. Only attribute with signatures that made large contribution
#    in step 1, above

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

# 3. Only use signatures previously observed in Lung-SCC

sigs.prop <- ExposureProportions("SBS96", "Lung-SCC")
extra <- c(SBS18 = 1/350, SBS29 = 3/350, SBS40 = 16/350, SBS44 = 1/350) # From vignette
sigs.prop <- c(sigs.prop, extra)

better1 <-
  apply(X          = lung.subset,
        MARGIN     = 2,
        FUN        = mSigAct::OptimizeExposureQP,
        signatures = s96[ , names(sigs.prop)])

better2 <- ICAMSxtra::SortExposure(better1)



ICAMSxtra::PlotExposureToPdf(better2, file = "better.ex.labels.pdf")

colnames(better2) <- NULL

ICAMSxtra::PlotExposureToPdf(better2, file = "better.ex.pdf")


# 4. Use the maximum a posteriori approach

mm <- mSigAct::DefaultManyOpts()
mm$trace <- 100
mm$global.opts$maxeval <- 10000

res <- tibble(sig.id = names(sigs.prop))
other.res <- list()

for (sampname in colnames(lung.subset)) {
  mapout <- MAPAssignActivity1(spect = lung.subset[ , sampname, drop = F],
                               sigs = s96,
                               sigs.presence.prop = sigs.prop,
                               max.level = 100,
                               p.thres = 0.01,
                               eval_f                  = ObjFnBinomMaxLHRound,
                               eval_g_ineq             = g_ineq_for_ObjFnBinomMaxLH2,
                               max.mc.cores = 50,
                               m.opts = mm
  )
  new.col <- mapout$MAP
  colnames(new.col)[2] <- sampname
  res <- dplyr::full_join(res, new.col)
  other.res <- c(other.res, sampname = mapout) # this is wrong, get literal sampname; next time use list(sampname, mapout)

}

save(other.res, file = "saved.other.res.Rdata")
# save(res, file = "saved.res.Rdata")

mr <- as.matrix(res[ , 2:ncol(res)])
rownames(mr) <- res$sig.id


all(is.na(mr["SBS44", ]))
all(is.na(mr["SBS29", ]))
all(is.na(mr["SBS18", ]))

all(is.na(mr["SBS40", ]))
all(is.na(mr["SBS8", ]))

xr <- mr[-c(10, 8, 7), ]

sum(xr, na.rm = T)
sum(mr, na.rm = T)

foo <- xr
foo[which(is.na(xr))] <- 0

foo2 <- ICAMSxtra::SortExposure(foo)

ICAMSxtra::PlotExposureToPdf(foo2, file = "MAP.ex.labels.pdf")

colnames(foo) <- NULL

ICAMSxtra::PlotExposureToPdf(foo2, file = "MAP.ex.pdf")

