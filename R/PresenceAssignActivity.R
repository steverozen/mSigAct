#' Find minimal set of signatures that can explain multiple spectra by first using
#' signature presence test
#'
#' @inheritParams MAPAssignActivity
#' 
#' @param sig.pres.test.nbinom.size The dispersion parameter for the negative
#'   binomial distribution used when conducting signature presence test first to
#'   filter out those signatures that are not needed in the reconstruction of
#'   the spectrum; smaller is more dispersed. See
#'   \code{\link[stats]{NegBinomial}}. If \code{NULL}, then use multinomial
#'   likelihood when conducting signature presence test.
#'   
#' @param sig.pres.test.p.thresh If the p-value of signature presence test for
#'   one signature is >= \code{sig.pres.test.p.thresh}, then this signature will
#'   not be used for assignment later.
#'   
#' @return A list with the elements:
#'
#' * \code{proposed.assignment}: The proposed set of signatures that can
#' plausibly explain \code{spectra}.
#' 
#' * \code{proposed.reconstruction}: The reconstruction based on
#' \code{proposed.assignment}.
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
#'  \code{PresenceAssignActivity} for each sample in \code{spectra}.
#' 
#' * \code{error.messages}: Only appearing if there are errors running
#' \code{PresenceAssignActivity}.
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
#' sigs.to.use <- SBS96.sigs[, names(sigs.prop), drop = FALSE]
#' retval <- PresenceAssignActivity(spectra = spectra,
#'                                  sigs = sigs.to.use,
#'                                  output.dir = file.path(tempdir(), "Lung-AdenoCA"),
#'                                  max.level = ncol(sigs.to.use) - 1,
#'                                  p.thresh = 0.05 / ncol(sigs.to.use),
#'                                  num.parallel.samples = 2,
#'                                  mc.cores.per.sample = 30,
#'                                  seed = 2561)
#'}
PresenceAssignActivity <- 
  function(spectra,
           sigs,
           output.dir,
           max.level                 = 5,
           p.thresh                  = 0.05,
           m.opts                    = DefaultManyOpts(),
           num.parallel.samples      = 5,
           mc.cores.per.sample       = min(20, 2^max.level),
           progress.monitor          = NULL,
           seed                      = NULL,
           drop.low.mut.samples      = TRUE,
           use.forward.search        = FALSE,
           save.files                = TRUE) {
    
    retval <- 
      MAPAssignActivity(spectra                   = spectra,
                        sigs                      = sigs,
                        use.sparse.assign         = TRUE,
                        use.forward.search        = use.forward.search,
                        output.dir                = output.dir,
                        max.level                 = max.level,
                        p.thresh                  = p.thresh,
                        m.opts                    = m.opts,
                        num.parallel.samples      = num.parallel.samples,
                        mc.cores.per.sample       = mc.cores.per.sample,
                        progress.monitor          = progress.monitor,
                        seed                      = seed,
                        max.subsets               = .Machine$double.xmax,
                        drop.low.mut.samples      = drop.low.mut.samples,
                        use.sig.presence.test     = TRUE,
                        save.files                = save.files)
    return(retval)
  }