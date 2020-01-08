cossim <- function(v1, v2) {
  if (!is.null(ncol(v1)))  {
    stopifnot(ncol(v1) == 1)
    v1 <- v1[ , 1]
  }
  if (!is.null(ncol(v2)))  {
    stopifnot(ncol(v2) == 1)
    v2 <- v2[ , 1]
  }
  lsa::cosine(v1, v2)
  
}