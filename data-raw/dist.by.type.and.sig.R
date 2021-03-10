distxx <- function(sig.name, cancer.type) {
  xx.by.type <- SplitPCAWGMatrixByTumorType(exposure$PCAWG$SBS96)
  ff <- xx.by.type[[cancer.type]]
  prop.0 <- sum(ff[sig.name, ] < 1) / ncol(ff)
  counts.non.0 <- ff[sig.name, ff[sig.name, ] >= 1]
  # par(mfrow = c(2, 1))
  if (length(counts.non.0) >= 10) {
    hist(counts.non.0, 
         breaks = seq(from = 0, to = max(counts.non.0) + 500, by = 500),
         main = paste(sig.name, cancer.type))
    plot(density(counts.non.0))
  }
  return(c(prop.0 = prop.0))
}

bb <- function() {
  pdf("foo.pdf")
  par(mfcol = c(6, 1))
  for (sig in rownames(PCAWG7::exposure$PCAWG$SBS96)) {
    for (cancer.type in CancerTypes()) {
      
      distxx(sig, cancer.type)
      
    }
    
  }
  dev.off()
}
