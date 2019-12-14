#' Load spectra rceived from Arnoud on 2019 Oct 23
#' 

LoadArnoudMCF10HepG2 <- function() {
  env <- new.env()
  load(devtools::package_file(
    file.path(
    "data-raw",
    "spectra.for.background.signatures",
    "MCF-10A-HepG2-background",
    "background_spectra.Rdata")),
    envir = env,
    verbose = TRUE)
  HepG2.background.spectra <- env$catSBS$catSBS96[ , 1:3]
  Mcf10A.background.spectra <- env$catSBS$catSBS96[ , 4:6]
  
  
  
}


#' Build a signature for background extraction from a matrix of spectra.
#' 
#' This function not only produces a signature, but also an
#' estimate of the number of mutations usually generated
#' by the signature and an indication of variability around
#' that estimate.
#' 
#' Only works on SBS 96 signatures.
#' 
#' @param spectra An \code{\link[ICAMS]{ICAMS}} \code{catalog} with 
#' \code{catalog.type = "counts"}.
#' 
#' @param algorithm See \code{\link[nloptr]{nloptr}}.
#' @param maxeval See \code{\link[nloptr]{nloptr}}.
#' @param print_level See \code{\link[nloptr]{nloptr}}.
#' @param xtol_rel See \code{\link[nloptr]{nloptr}}.
#' @param xtol_abs See \code{\link[nloptr]{nloptr}}.
#' 
#' @return A list with the elements 
#' \enumerate{
#' \item \code{signature} An \code{\link[ICAMS]{ICAMS}}
#' \code{catalog} with
#' \code{catalog.type == "counts.signature"}.
#' \item \code{log10.counts} Mean log base 10 of the 
#' total counts in \code{spectra}
#' \item \code{sd.log10.counts.per.base} Standard deviation of 
#' \code{log10.counts.per.base}.
#' }
#' 
#' @export
#
# Internal notes: this function could estimate the
# dispersion for each channel separately, or we could
# even find the maximum likelihood estimate of the
# generating signature, the number of 
# counts it produces, and its dispersion parameter.

EstimateSignatureFromSpectraLH <-
  function(spectra,
           algorithm="NLOPT_LN_COBYLA",
           maxeval=1000, 
           print_level=0,
           xtol_rel=0.001,  # 0.0001,)
           xtol_abs=0.0001)
  {
    # Maybe this is excessively "realistic"; maybe
    # just take mean of spectra.
    # It turns out that the result signature is exactly the mean.
    # 
    # But this can estimate the negative binomial dispersion
    # parameter too, so perhaps useful.
    
    spectra.as.sigs <-
      ICAMS::TransformCatalog(spectra, 
                              target.catalog.type = "counts.signature")
    x0.sig.vec <- rowSums(spectra.as.sigs) / ncol(spectra)
    # Start with a signature that is an average
    
    mean.sig <- matrix(x0.sig.vec, ncol = 1)
    rownames(mean.sig) <- rownames(spectra)
    mean.sig <- ICAMS::as.catalog(
      mean.sig,
      ref.genome   = attr(spectra, "ref.genome",   exact = TRUE),
      region       = attr(spectra, "region",       exact = TRUE),
      catalog.type = "counts.signature")
    
    x0.sig.and.size <- c(x0.sig.vec, 200)
    
    
    ret <- nloptr::nloptr(
      x0 = x0.sig.and.size,
      eval_f = NegLLHOfSignature,
      lb = rep(0, length(x0.sig.and.size)),
      opts = list(algorithm   = algorithm,
                  xtol_rel    = xtol_rel,
                  print_level = print_level,
                  maxeval     = maxeval),
      spectra = spectra)
    
    len <- nrow(spectra)
    sig <- matrix(ret$solution[1:len], ncol = 1)
    rownames(sig) <- ICAMS::catalog.row.order$SBS96
    sig <- sig / sum(sig)
    sig <- ICAMS::as.catalog(
      sig,
      region       = attr(spectra.as.sigs, "region",     exact = TRUE),
      ref.genome   = attr(spectra.as.sigs, "ref.genome", exact = TRUE),
      catalog.type = "counts.signature")
    
    nbinom.size <- ret$solution[len + 1]
    
    return(list(background.sig   = sig,
                mean.sig         = mean.sig,
                nbinom.size      = nbinom.size,
                count.nbinom.mu  = mean(colSums(spectra)),
                count.binom.size = 20,
                nloptr.ret       = ret))
    
  }



# Part 1, Background spectra for HepG2.

#' Make a standard \code{\link[ICAMS]{ICAMS}} SBS 96 catalog for the HepG2 background signature.
#' 
#' @keywords internal

MakeHepG2BackgroundPart1 <- function() {
  spectra <-
    ICAMS:::ReadDukeNUSCat192(
      file = system.file(
        "data-raw/spectra.for.background.signatures/HepG2_SC2_background.txt",
        package = "mSigAct"),
      catalog.type = "counts",
      region = "unknown")$cat96
  ICAMS::WriteCatalog(
    spectra,
    file = system.file(
      "data-raw/spectra.for.background.signatures/HepG2.background.96.csv",
      package = "mSigAct")
  )
}


#' Read the specified spectra file and estimate the HepG2 background signature.
#' 
#' @keywords internal
#' 
MakeHepG2BackgroundPart2 <- function(maxeval) {
  set.seed(3214)
  s1 <- ICAMS::ReadCatalog(
    system.file(
      "data-raw/spectra.for.background.signatures/HepG2.background.96.csv",
      package = "mSigAct"),
    region = "genome",
    catalog.type = "counts")
  
  # Calculate parameters for a negative bionomial 
  # distribution modeling the
  # number of mutations generated by the signature.
  count.nbinom.mu <- mean(colSums(s1))
  count.nbinom.size <- 20 # determined manually TODO(Steve): analytical formula?
  
  
  ret2 <-
    EstimateSignatureFromSpectraLH(s1, maxeval = maxeval, print_level = 1)
  
  return(list(background.sig    = ret2[["background.sig"]],
              mean.sig          = ret2[["mean.sig"]],
              sig.nbinom.size   = ret2[["nbinom.size"]],
              count.nbinom.mu   = count.nbinom.mu, 
              count.nbinom.size = count.nbinom.size))
}


# Part 2, Spectra from cisplatin exposed HepG2 cells.

#' Make spectrum catalog from VCFs from cisplatin exposed HepG2
#' @keywords internal
#' @return The catalog
MakeCisplatinCatalogs <- function() {
  files <- dir("tests/testthat/test.data/HepG2_Cis/", full.names = TRUE)
  cats <- ICAMS::StrelkaSBSVCFFilesToCatalog(
    files = files,
    ref.genome = "hg19",
    trans.ranges = 
      ICAMS::trans.ranges.GRCh37,
    region = "genome")
  ICAMS::WriteCatalog(cats$catSBS96, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS96.csv")
  ICAMS::WriteCatalog(cats$catSBS192, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS192.csv")
  ICAMS::WriteCatalog(cats$catDBS78, 
                      file = "tests/testthat/test.data/HepG2_Cis/DBS78.csv")
  ICAMS::PlotCatalogToPdf(cats$catSBS96, 
                          file = "tests/testthat/test.data/HepG2_Cis/SBS96.pdf")
  ICAMS::PlotCatalogToPdf(cats$catSBS192, 
                          file = "tests/testthat/test.data/HepG2_Cis/SBS192.pdf")
  ICAMS::PlotCatalogToPdf(cats$catDBS78, 
                          file = "tests/testthat/test.data/HepG2_Cis/DBS78.pdf")
  
  return(cats)
}

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

# mcf10a <- LoadToEnvironment(
#   "data-raw/spectra.for.background.signature/MCF-10A-background/background_spectra.Rdata")
