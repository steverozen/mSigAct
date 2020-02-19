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
