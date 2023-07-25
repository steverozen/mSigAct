#' Find minimal set of signatures that can explain multiple spectra by first using
#' signature presence test
#'
#' @inheritParams MAPAssignActivity
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
#' * \code{time.for.assignment}: Value from \code{system.time} for running
#'  \code{PresenceAssignActivity} for each sample in \code{spectra}.
#' 
#' * \code{error.messages}: Error messages running
#' \code{PresenceAssignActivity}.
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
#'                                  num.parallel.samples = 2,
#'                                  mc.cores.per.sample = 30,
#'                                  seed = 2561)
#'}
PresenceAssignActivity <- 
  function(spectra,
           sigs,
           output.dir,
           p.thresh                  = DefaultPThresh(sigs),
           m.opts                    = DefaultManyOpts(spectra = spectra),
           num.parallel.samples      = 1,
           mc.cores.per.sample       = 1,
           seed                      = 123,
           drop.low.mut.samples      = FALSE,
           save.files                = TRUE) {
    
    retval <- 
      MAPAssignActivity(spectra                   = spectra,
                        sigs                      = sigs,
                        use.forward.search        = TRUE,
                        output.dir                = output.dir,
                        p.thresh                  = p.thresh,
                        m.opts                    = m.opts,
                        num.parallel.samples      = num.parallel.samples,
                        mc.cores.per.sample       = mc.cores.per.sample,
                        seed                      = seed,
                        drop.low.mut.samples      = drop.low.mut.samples,
                        use.sig.presence.test     = TRUE,
                        save.files                = save.files)
    return(retval)
  }