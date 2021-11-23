#' Calculate various distances information between spectra and reconstructed spectra
#'
#' @inheritParams AddSigActivity
#'
#' @return A data frame with row names being the sample names in \code{spectra}.
#' The columns show various distances information between \code{spectra} and reconstructed
#' spectra using \code{sigs} and \code{exposure}. The column names are:
#' * log.likelihood: `log(likelihood(spectra | reconstructed spectra))`
#' * euclidean: Euclidean distance
#' * manhattan: Manhattan distance
#' * cosine: Cosine similarity
#' * mutations: Mutations of `spectra`
#' * scaled.manhattan: `manhattan` / `mutations`
#' * scaled.euclidean: `euclidean` / `mutations`
#' 
#' @md
#' 
#' @keywords internal
CalculateDistance <- 
  function(spectra, exposure, sigs, sigs.presence.prop, nbinom.size = 5,
           likelihood.dist = "multinom", use.sparse.assign = TRUE) {
    # Check whether there are some samples which have zero mutations
    retval <- RemoveZeroMutationSample(spectra = spectra, exposure = exposure)
    spectra <- retval[["spectra"]]
    exposure <- retval[["exposure"]]
    
    sig.activities <- AddSigActivity(spectra = spectra,
                                     exposure = exposure,
                                     sigs = sigs,
                                     sigs.presence.prop = sigs.presence.prop,
                                     nbinom.size = nbinom.size,
                                     likelihood.dist = likelihood.dist,
                                     use.sparse.assign = use.sparse.assign)
    distance.list <- lapply(sig.activities, FUN = function(x) {
      return(x$distances[, 2, drop = FALSE])
    })
    distances <- do.call("cbind", distance.list)
    colnames(distances) <- colnames(spectra)
    rownames(distances) <- sig.activities[[1]]$distances$method
    
    tmp <- rbind(distances, colSums(spectra))
    rownames(tmp)[5] <- "mutations"
    
    distance.info <- as.data.frame(t(tmp))
    distance.info$scaled.manhattan <- 
      distance.info$manhattan / distance.info$mutations
    distance.info$scaled.euclidean <- 
      distance.info$euclidean / distance.info$mutations
    return(distance.info)
  }