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
    return(RemoveZeroActivitySig(list.of.exposures[[1]]))
  }
  
  tmp <- lapply(list.of.exposures, FUN = function(one.exposure) {
    return(as.data.frame(t(one.exposure)))
  })
  
  combined.exposure <- t(dplyr::bind_rows(tmp))
  combined.exposure[is.na(combined.exposure)] <- 0
  combined.exposure <- 
    combined.exposure[SortSigId(rownames(combined.exposure)), , drop = FALSE]
  
  return(RemoveZeroActivitySig(combined.exposure))
}