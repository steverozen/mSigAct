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