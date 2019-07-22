#' @keywords internal
LikelihoodBackgroundOnly <- function(spectra, background.sig) {
     # Find max likelihood and count to generate each
     # spectrum in spectra using only background (the model could contain 
     # an expected distribution of intensity of background mutagenesis)
  
     # will compare this with a model with another signature
     # (the same across all samples) with variable count
     # per sample
  
     # To find the background signature, find the shape and
     # and intensity that gives the observed backgrounds.
  
}



#' Subtract background determined by maxium likelihood of
#' spectrum1 = b1 * background + (total-mut1 - b1) * sig * prob(b1)
#' spectrum2 = b2 * background + (total-mut2 - b2) * sig * prob(b2)
#' ...
#' 
#' (Arnoud's idea - take low regions in spectrum and
#' assume they are due to background)
#' 
#' @param spectra An \code{ICAMS}
#'  \code{counts} or \code{density} \code{catalog}.
#'
#' @param background.sig An \code{ICAMS} 
#'  \code{counts.signature} or \code{density.signature} \code{catalog}. 
#'  
#'  \itemize{
#'  
#'  \item If
#'  the class and \code{background.sig} and \code{catalog}
#'  must be the same, and
#'  
#'  \item if \code{catalog.type} of
#'  \code{background.sig} is \code{density.signature} then the \code{catalog.type}
#'  of \code{spectra} must be \code{density}, and
#'  
#'  \item If the \code{catalog.type} of
#'  \code{background.sig} is \code{counts.signature} then the \code{catalog.type}
#'  of \code{spectra} must be \code{counts} \strong{and} the
#'  \code{abundance} attributes must match.
#'  
#'  } 
#'  
#' @export
LSubtractBackground <- function(spectra, background.sig) {
  out.spectra <- 
    apply(X = colnames(spectra),
          MARGIN = 2,
          function(spectrum) {
            SubtractBackground1(spectrum, background.sig, background.count)
          }
    )
  return(out.spectra) # Not sure how to make this a catalog again
}




#' Subtract background if required by a mutation count test methods
#' 
#' @inheritParams LSubtractBackground
#'
#' @param background.count Expected total number of mutatations due to 
#' \code{background.sig}.
#'  
#' @export
CSubtractBackground <- function(spectra, background.sig, background.count) {
    out.spectra <- 
      apply(X = colnames(spectra),
            MARGIN = 2,
            function(spectrum) {
              SubtractBackground1(spectrum, background.sig, background.count)
            }
      )
    return(out.spectra) # Not sure how to make this a catalog again
}



#' Estimate a signature from several spectra; two possible approaches,
#' 1. Just an average of the signatures corresponding to the spectra.
#' 2. Find a maximum likelihood signature and coefficients (exposures)
#' in the clones.
#  For now, we assume that option 1 is sufficient.

GetHepG2Background1 <- function() {
  # This is a 192-channel catalog, but it is whole genome,
  # so we use the collapsed version.
  c1 <- ReadDukeNUSCat192(
    devtools::package_file(
      "data-raw/possible-future-code/HepG2_SC2_clones.txt"),
    ref.genome = "hg19", region = "unknown",
    catalog.type = "counts")
  
  counts96 <- as.catalog(c1$cat96, ref.genome = "hg19",
                         region = "genome")
  
  PlotCatalogToPdf(
    counts96, 
    devtools::package_file("data-raw/possible-future-code/counts96.pdf"))
  PlotCatalogToPdf(
    TransformCatalog(counts96, target.catalog.type = "density"),
    devtools::package_file("data-raw/possible-future-code/density96.pdf"))
  
  pre.sig96 <-
    TransformCatalog(c1$cat96, target.catalog.type = "counts.signature")
  
  sig96 <-
    as.catalog(
      as.matrix(apply(pre.sig96, MARGIN = 1, mean), nrow = 96),
      ref.genome = "hg19",
      region = "genome",
      catalog.type = "counts.signature")
  
  PlotCatalogToPdf(
    sig96,
    devtools::package_file("data-raw/possible-future-code/sig96.pdf"))
  
  PlotCatalogToPdf(
    TransformCatalog(sig96, target.catalog.type = "density.signature"),
    devtools::package_file("data-raw/possible-future-code/sig96-dens.pdf"))
  
  return(sig96)
  
}

# https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/

