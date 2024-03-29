#' Find known signatures that can most sparsely reconstruct each spectrum in a catalog.
#'
#' @inheritParams MAPAssignActivity
#'
#' @return A list with the elements:
#'
#' * \code{proposed.assignment}: The most sparse set of signatures that can
#' plausibly explain \code{spectra}.
#' 
#' * \code{proposed.reconstruction}: The reconstruction based on sparse
#' assignment.
#'
#' * \code{reconstruction.distances}: Various distances and similarities
#' between \code{spectra} and \code{proposed.reconstruction}.
#' 
#' * \code{all.tested}: All tested possible ways to reconstruct each
#' sample in \code{spectra}.
#' 
#' * \code{alt.solutions}: A \code{tibble} showing all the alternative solutions
#' that are statistically as good as the \code{proposed.assignment} that can
#' plausibly reconstruct \code{spectra}.
#' 
#' * \code{time.for.assignment}: Value from \code{system.time} for running
#'  \code{SparseAssignActivity} for each sample in \code{spectra}.
#' 
#' * \code{error.messages}: Only appearing if there are errors running
#' \code{SparseAssignActivity}.
#'
#' The elements \code{proposed.assignment}, \code{proposed.reconstruction},
#' \code{reconstruction.distances}, \code{all.tested},
#' \code{time.for.assignment} will be \code{NULL} if the algorithm could not
#' find the optimal reconstruction or there are errors coming out for
#' \strong{all} samples.
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
#' SBS96.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
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
           max.level               = 5,
           p.thresh                = 0.05,
           m.opts                  = DefaultManyOpts(),
           num.parallel.samples    = 5,
           mc.cores.per.sample     = min(20, 2^max.level),
           progress.monitor        = NULL,
           seed                    = NULL,
           max.subsets             = 1000,
           drop.low.mut.samples    = TRUE) {
    
    retval <- MAPAssignActivity(spectra                 = spectra,
                                sigs                    = sigs,
                                use.sparse.assign       = TRUE,
                                output.dir              = output.dir,
                                max.level               = max.level,
                                p.thresh                = p.thresh,
                                m.opts                  = m.opts,
                                num.parallel.samples    = num.parallel.samples,
                                mc.cores.per.sample     = mc.cores.per.sample,
                                progress.monitor        = progress.monitor,
                                seed                    = seed,
                                max.subsets             = max.subsets,
                                drop.low.mut.samples    = drop.low.mut.samples)
    return(retval)
  }

