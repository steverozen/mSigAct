library(mSigAct)

source("data-raw/errors/ID.attribution.error/code/debug-ID-attribution.R")


MAP.out <- DebugIDAttribution()

plus.minus.grid <- function(spectrum, sigs) {
  lapply(rownames(spectrum), 
        function(nn) { c(spectrum[nn, ] / sigs[nn, ]) }
        
        )
  
  
}

file.ID <- "data-raw/errors/ID.attribution.error/data/catID.counts.csv"
catID <- ICAMS::ReadCatalog(file = file.ID, ref.genome = "hg19", region = "genome")
mm <- rmultinom(1000, size = sum(catID), prob = catID)
View(mm)
# each row is a channel, each column is a sample
colSums(mm)
sum(catID)
cos <- apply(mm, MARGIN = 2, FUN = function(mmm) cossim(mmm, catID))
dmultinom(catID[ , 1], prob = MAP.out[[1]]$proposed.reconstruction[ , 1], log = T)
mSigAct:::LLHSpectrumNegBinom(catID[ , 1],
                              MAP.out[[1]]$proposed.reconstruction[ , 1])
