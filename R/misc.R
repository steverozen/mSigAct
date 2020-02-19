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

ClosestCosSig <- function(spectrum) {
  cos <- 
    apply(PCAWG7::signature$genome$SBS96, 
          MARGIN = 2, 
          FUN = 
            function(sig) {
              lsa::cosine(as.vector(sig), as.vector(spectrum))})
  max.cos <- which(cos == max(cos))
  return(cos[max.cos])
  
}


ClosestCosSigDensity <- function(spectrum) {
  spec <- 
    ICAMS::TransformCatalog(
      spectrum,
      target.catalog.type = "density")
  
  sigs <- 
    ICAMS::TransformCatalog(
      PCAWG7::signature$genome$SBS96,
      target.catalog.type = "density.signature")
  
  cos <- 
    apply(sigs, 
          MARGIN = 2, 
          FUN = 
            function(sig) {
              lsa::cosine(as.vector(sig), as.vector(spec))})
  max.cos <- which(cos == max(cos))
  return(cos[max.cos])
  
}


LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

ClearWarnings <- function() {
  assign("last.warning", NULL, envir = baseenv())
  
}
