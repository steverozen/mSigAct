#' Generate synthetic data based on control data from Kucab et al., 2019.
#'
#' @param target.sig The sig to test for, in Kucab's format.
#' @param target.ratio.to.control The number of target signature
#'    mutations as a proporiton of the number of background signatures.
#'    E.g. 1 means the same number of mutations due to the target signatures
#'    as due to the background signature; 2 means twice as many mutations
#'    due to the the target signature as due to the background signature.
#' @param num.replicates Number of synthetic signatures to generate.
#' @param seed Optional random seed.
#' @param add.noise If \code{TRUE} and negative binomial noise to the
#'    proportion of the spectra due to the target signature.
#'    
#'  The background contributions are taken randomly from the control
#'  spectra in Kucab et al., 2019.
Generate1KucabSynData <- 
  function(target.sig, 
           target.ratio.to.control, 
           num.replicates, 
           seed = NULL,
           add.noise = TRUE) {
  if (is.null(seed)) seed <- 10101
  set.seed(seed,"L'Ecuyer-CMRG")
  controls <- sample(35, size = num.replicates, replace = TRUE)
  control.spectra <- 
    mSigAct::kucab.control.spectra.kucab[ , controls, drop = FALSE]
  stopifnot(rownames(target.sig) == rownames(control.spectra))

  num.control.muts <- colSums(control.spectra)
  
  num.target.muts <- num.control.muts * target.ratio.to.control
  
  target.spectra.no.noise <- 
    as.matrix(target.sig) %*% matrix(num.target.muts, nrow = 1)
  
  if (add.noise) {
    target.spectra.noise <- 
      matrix(sapply(
        round(target.spectra.no.noise),
        function(mu) stats::rnbinom(n = 1, size = 10000, mu = mu)),
        nrow = nrow(target.spectra.no.noise))
    
    test.spectra <- control.spectra + target.spectra.noise
  } else {
    test.spectra <- control.spectra + target.spectra.no.noise
  }
  colnames(test.spectra) <- paste0("test.spec.", 1:ncol(test.spectra))
  test.spectra <- round(test.spectra)

  # Need to return the exposures too
  return(list(test.spectra = test.spectra,
              bg.exposures = num.control.muts,
              target.exposures = num.target.muts))
}
  
#' Convert a Kucab et al., spectra catalog to an ICAMS catalog.
#' 
#' The Zou code in Kucab et al., 2019 uses an unusual row order
#' and the first column has the mutation names as e.g.
#' A[C>T]G.
#
#' @param m The spectra (or signature) as a matrix, with
#'  rownames e.g. "A[C>A]A", "A[C>G]C", ...
#'  
#' @export
#'  
KucabToICAMSSpectra <- function(m) {
  stopifnot(rownames(m) == rownames(mSigAct::kucab.sub.catalog))
  stopifnot(is.numeric(m[ ,1]))
  new.rownames <- ICAMS:::Unstaple96(rownames(m))
  rownames(m) <- new.rownames
  m <- m[ICAMS::catalog.row.order$SBS96, ]
  return(ICAMS::as.catalog(m, region = "genome", catalog.type = "counts"))
}
