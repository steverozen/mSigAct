#' Generate exposure model priors from exposure data for one cancer type.
#'
#' @param exposure.matrix A matrix or numeric data frame in which each row is a signature
#'  (with signature identifiers as rownames) and each column is a 
#'  sample (with colnames as sample identifiers)
#' 
GenerateExposureModelOneTumorType <- function(exposure.matrix, cancer.type) {
  for (sig.name in rownames(exposure.matrix)) {
    exposures <- exposure.matrix[sig.name, ]
    exp2 <- exposures[exposures > 0]
    if (length(exp2) == 0) {
      cat("No exposure to ", sig.name, "\n")
      next
      
    }
    # old.digits <- getOption("digits")
    # old.scipen <- getOption("scipen")
    # on.exit(options(digits = old.digits, scipen =  old.scipen))
    # options(digits = 3, scipen = 2)
    prop.non.zero <- length(exp2) / length(exposures)
    hist(exp2, 
         main = paste(
           cancer.type, 
           sig.name, "num non-0 = ", length(exp2), ", proportion non-0 =", format(prop.non.zero, digits = 3)), 
         breaks = seq(from = 0, to = max(exp2) + 100, by = 100),
         xlab = "Number of mutations")
    
  }

}



SplitPCAWG7ExposureByType <- function(exp) {
  types <- strsplit(colnames(exp), "::", fixed = TRUE)
  types <- unlist(lapply(types, "[", 1))
  exp2 <- split(as.data.frame(t(exp)), types)  
  
  exp3 <- lapply(exp2, t)
  
  return(exp3)
  
}

if (FALSE) {
  exp.by.type <- SplitPCAWG7ExposureByType(PCAWG7::exposure$PCAWG$SBS96)
  
  pdf("data-raw/PCAWG.exposures.by.sig.pdf")
  par(mfrow = c(4, 1))
  
  for (cancer.type in names(exp.by.type)) {
    cat(cancer.type, "\n")
    foo <- GenerateExposureModelOneTumorType(exp.by.type[[cancer.type]], cancer.type)
    
  }
  
  dev.off()

  
}