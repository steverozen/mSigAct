DefaultPThresh <- function(sigs) {
  p.thresh <- 0.001 / (4 * ncol(sigs))
  return(p.thresh)
}