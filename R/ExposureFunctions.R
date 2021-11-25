#' @keywords internal
MergeTwoExposures <- function(exposure1, exposure2) {
  if (inherits(exposure1, "numeric")) {
    exposure1 <- as.data.frame(exposure1)
  }
  
  if (inherits(exposure2, "numeric")) {
    exposure2 <- as.data.frame(exposure2)
  }
  tmp <- dplyr::bind_rows(as.data.frame(t(exposure1)), 
                          as.data.frame(t(exposure2)))
  tmp[is.na(tmp)] <- 0
  
  merged.exp <- t(tmp)
  merged.exp2 <- merged.exp[SortSigId(rownames(merged.exp)), , drop = FALSE]
  return(merged.exp2)
}

#' @keywords internal
MergeListOfExposures <- function(list.of.exposures) {
  num.of.exposure <- length(list.of.exposures)
  
  if (num.of.exposure < 2) {
    stop("Only one exposure, no need to merge")
  }
  
  combined.exposure <- 
    MergeTwoExposures(exposure1 = list.of.exposures[[1]],
                      exposure2 = list.of.exposures[[2]])
  
  if (num.of.exposure > 2) {
    for (i in 3:num.of.exposure) {
      combined.exposure <- 
        MergeTwoExposures(exposure1 = combined.exposure,
                          exposure2 = list.of.exposures[[i]])
    }
  } 
  return(RemoveZeroActivitySig(combined.exposure))
}