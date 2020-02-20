# To create background information variables, at top level, run
# the following functions
# 
# MakeMCF10HepG2BackgroundVars()
# SaveHepG2andMCF10ABackgroundInfo()

#' Create spectra from VCFs and load spectra as package variables.
#' VCFs received from Arnoud on 2020 Feb 04


MakeMCF10HepG2BackgroundVars <- function() {
  data.dir <- file.path("data-raw",
                        "background.signature.spectra", 
                        "backgroundVCFS-2020-02-04")
  sbs.files <- dir(path = data.dir, pattern = "SNV", full.names = TRUE)
  names <- sub(".*((MCF10A|HepG2)_SC._cl.).*", "\\1", sbs.files, perl = TRUE)
  
  # No substantial indel background, so will not read these VCFs
  # id.files  <- dir(path = data.dir, pattern = "INDEL", full.names = TRUE)
  sbs.cat <- 
    ICAMS::StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
      files         = sbs.files, 
      ref.genome    = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
      region        = "genome",
      names.of.VCFs = names,
      output.file = file.path(data.dir, "MCF10A-and-HepG2-bg"))
  ICAMS::WriteCatalog(sbs.cat$catSBS96, 
                      file.path(data.dir, "MCF10A-and-HepG2-bg-SBS96-bg.csv"))
  ICAMS::WriteCatalog(sbs.cat$catSBS192, 
                      file.path(data.dir, "MCF10A-and-HepG2-bg-SBS192-bg.csv"))
  
  if (FALSE) {
    # Old code
    env <- new.env()
    load(devtools::package_file(
      file.path(
        "data-raw",
        "spectra.for.background.signatures",
        "MCF-10A-HepG2-background",
        "background_spectra.Rdata")),
      envir = env,
      verbose = TRUE)
  }
  
  # The ref genomes have an arbitrary file path encoded in them
  attr(sbs.cat$catSBS96, "ref.genome") <- NULL

  HepG2.background.spectra  <- sbs.cat$catSBS96[ , 1:3]
  MCF10A.background.spectra <- sbs.cat$catSBS96[ , 4:6]
  usethis::use_data(HepG2.background.spectra, overwrite = TRUE)
  usethis::use_data(MCF10A.background.spectra, overwrite = TRUE)
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
                              target.catalog.type = "counts.signature",
                              target.abundance    = 
                                attr(spectra, "abundance", exact = TRUE))
    x0.sig.vec <- rowSums(spectra.as.sigs) / ncol(spectra)
    # Start with a signature that is an average
    
    mean.sig <- matrix(x0.sig.vec, ncol = 1)
    rownames(mean.sig) <- rownames(spectra)
    mean.sig <- ICAMS::as.catalog(
      mean.sig,
      ref.genome   = NULL,
      abundance    = attr(spectra, "abundance",    exact = TRUE),
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
      abundance    = attr(spectra.as.sigs, "abundance",  exact = TRUE),
      ref.genome   = NULL,
      catalog.type = "counts.signature")
    
    nbinom.size <- ret$solution[len + 1]
    
    return(list(background.sig   = sig,
                mean.sig         = mean.sig,
                nbinom.size      = nbinom.size,
                count.nbinom.mu  = mean(colSums(spectra)),
                count.binom.size = 20,
                nloptr.ret       = ret))
    
  }


#' Estimate a background signature.
#' 
#' @keywords internal
#' 
MakeBackground <- function(bg.spectra, maxeval) {
  set.seed(3214)

  # Calculate parameters for a negative bionomial 
  # distribution modeling the
  # number of mutations generated by the signature.
  count.nbinom.mu <- mean(colSums(bg.spectra))
  count.nbinom.size <- 20 # determined manually TODO(Steve): analytical formula?
  
  
  ret2 <-
    EstimateSignatureFromSpectraLH(
      bg.spectra,
      maxeval = maxeval, 
      print_level = 1)
  
  return(list(background.sig    = ret2[["background.sig"]],
              mean.sig          = ret2[["mean.sig"]],
              sig.nbinom.size   = ret2[["nbinom.size"]],
              count.nbinom.mu   = count.nbinom.mu, 
              count.nbinom.size = count.nbinom.size))
}

SaveHepG2andMCF10ABackgroundInfo <- function() {
  HepG2.background.info <-
    MakeBackground(mSigAct::HepG2.background.spectra, 10000)
  MCF10A.background.info <-
    MakeBackground(mSigAct::MCF10A.background.spectra, 10000)
  usethis::use_data(HepG2.background.info, overwrite = TRUE)
  usethis::use_data(MCF10A.background.info, overwrite = TRUE)
}

