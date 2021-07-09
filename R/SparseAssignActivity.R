#' Find known signatures that can most sparsely reconstruct each spectrum in a catalog.
#'
#' @param spectra The spectra (multiple spectra) to be reconstructed.
#'
#' @param sigs The known signatures to use in reconstruction.
#' 
#' @param output.dir Directory path to save the output file.
#'
#' @param max.level The largest number of signatures to consider discarding
#' in the reconstruction.
#'
#' @param p.thresh The maximum p value based on which it is decided
#' to retain a signature in a reconstruction.
#'
#' @param m.opts For documentation
#'    see \code{\link{DefaultManyOpts}}.
#'
#' @param num.parallel.samples The (maximum) number of samples to run in parallel; each
#'    sample in turn can require multiple cores, as governed by
#'    \code{mc.cores.per.sample}.
#'
#' @param mc.cores.per.sample
#'    The maximum number of cores to use for each sample.
#'    On Microsoft Windows machines it is silently changed to 1.
#'    
#' @param seed Random seed; set this to get reproducible results. (The
#'   numerical optimization is in two phases; the first, global phase
#'   might rarely find different optima depending on the random
#'   seed.)
#'
#' @return A list with the elements:
#'
#' * \code{proposed.assignment}: Proposed sparse assignment for \code{spectra}.
#'
#' * \code{proposed.reconstruction}: Proposed reconstruction of \code{spectra}
#' based on sparse assignment.
#'
#' * \code{reconstruction.distances}: Various distances and similarities
#' between \code{spectra} and \code{proposed.reconstruction}.
#'
#' @md
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # This is a long running example unless parallel computing is supported on your machine
#' indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
#' spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
#' sigs <- PCAWG7::signature$genome$SBS96
#' sigs.prop <- ExposureProportions(mutation.type = "SBS96",
#'                                  cancer.type = "Lung-AdenoCA")
#' sigs.to.use <- sigs[, names(sigs.prop), drop = FALSE]
#' sparse.out <- SparseAssignActivity(spectra = spectra,
#'                                    sigs = sigs.to.use,
#'                                    output.dir = file.path(tempdir(), "Lung-AdenoCA"),
#'                                    max.level = ncol(sigs.to.use) - 1,
#'                                    p.thresh = 0.05 / ncol(spectra),
#'                                    num.parallel.samples = 2,
#'                                    mc.cores.per.sample = 30,
#'                                    seed = 2561)
#'}
SparseAssignActivity <-
  function(spectra,
           sigs,
           output.dir,
           max.level            = 5,
           p.thresh             = 0.05,
           m.opts               = DefaultManyOpts(),
           num.parallel.samples = 5,
           mc.cores.per.sample  = min(20, 2^max.level),
           seed                 = NULL) {
    
    
    f1 <- function(i) {
      retval1 <- RunSparseAssignOneSample(
        spect                   = spectra[ , i, drop = FALSE],
        sigs                    = sigs,
        output.dir              = output.dir,
        max.level               = max.level,
        p.thresh                = p.thresh,
        m.opts                  = m.opts,
        max.mc.cores            = mc.cores.per.sample,
        seed                    = seed)
      
      return(retval1)
    }
    
    num.parallel.samples <- Adj.mc.cores(num.parallel.samples)
    
    retval <- parallel::mclapply(1:ncol(spectra),
                                 f1,
                                 mc.cores = num.parallel.samples)
    
    names(retval) <- colnames(spectra)
    
    proposed.assignment <- GetExposureInfo(list.of.MAP.out = retval)
    # Replace NA to 0 in proposed.assignment
    proposed.assignment[is.na(proposed.assignment)] <- 0
    
    proposed.reconstruction <- GetReconstructionInfo(list.of.MAP.out = retval)
    # Add attributes to proposed.reconstruction to be same as spectra
    proposed.reconstruction <- AddAttributes(proposed.reconstruction, spectra)
    
    reconstruction.distances <- GetDistanceInfo(list.of.MAP.out = retval, 
                                                sparse.assign = TRUE)
    
      return(list(proposed.assignment          = proposed.assignment,
                  proposed.reconstruction      = proposed.reconstruction,
                  reconstruction.distances     = reconstruction.distances))
    
    if (FALSE) {
      f1 <- function(i) {
        retval1 <- SparseAssignActivity1(
          spect        = spectra[ , i, drop = FALSE],
          sigs         = sigs,
          p.thresh     = p.thresh,
          m.opts       = m.opts,
          max.level    = max.level,
          max.mc.cores = mc.cores.per.sample,
          seed         = seed)
        
        return(retval1)
      }
      
      if (is.null(m.opts)) m.opts <- DefaultManyOpts()
      
      num.parallel.samples <- Adj.mc.cores(num.parallel.samples)
      
      retval <- parallel::mclapply(1:ncol(spectra),
                                   f1,
                                   mc.cores = num.parallel.samples)
      check.mclapply.result(
        retval, "SparseAssignActivity", colnames(spectra))
      
      other.info <- lapply(retval, attributes)
      retval2 <- matrix(unlist(retval), ncol = length(retval))
      colnames(retval2) <- colnames(spectra)
      rownames(retval2) <- colnames(sigs)
      
      return(list(exposure = retval2, other.info = other.info))
    }
}

