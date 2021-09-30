#' Run \code{\link{MAPAssignActivity1}} on one sample from the PCAWG platinum data set.
#'
#' @param cancer.type A cancer type from the PCAWG exposures matrix.
#'
#' @param sample.index The index of the sample within the exposures matrix.
#'
#' @param mutation.type One of "SBS96", "SBS192", "ID", "DBS78"
#'
#' @param max.level The maximum number of signatures to try removing.
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'
#' @param out.dir If non-NULL create this directory if necessary and put
#'   results there.
#'
#' @param max.mc.cores
#'   The maximum number of cores to use.
#'   On Microsoft Windows machines it is silently changed to 1.
#'
#' @param p.thresh If
#'  the p value for a better reconstruction with than without a set of signatures
#'  is > than \code{p.thresh}, then we can use exposures without this set.
#'
#' @param max.presence.proportion  The maximum value of the proportion
#'   of tumors that must have a given signature. Used so that it is
#'   possible to exclude a signature from a spectrum, e.g.
#'   perhaps all examples of tumor types have SBS5, but we want
#'   to allow a small chance that SBS5 is not present.
#'
#' @param sigs.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being signature identifiers. Can be the
#'    return value from \code{\link{ExposureProportions}}.
#'
#' @return See \code{\link{OneMAPAssignTest}}.
#'
#' @keywords internal
#'

PCAWGMAPTest <- function(cancer.type,
                         sample.index,
                         mutation.type,
                         max.level = 5,
                         max.mc.cores,
                         m.opts = DefaultManyOpts(),
                         out.dir = NULL,
                         p.thresh = 0.01,
                         max.presence.proportion = 0.99,
                         sigs.prop               = NULL) {

  exposure.mutation.type <-
    ifelse(mutation.type == "SBS192", "SBS96", mutation.type)

  stopifnot(exposure.mutation.type %in% names(PCAWG7::exposure$PCAWG))

  p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(
    PCAWG7::spectra$PCAWG[[mutation.type]])
  stopifnot(cancer.type %in% names(p7))

  ex <- PCAWG7::SplitPCAWGMatrixByTumorType(
    PCAWG7::exposure$PCAWG[[exposure.mutation.type]])
  stopifnot(cancer.type %in% names(ex))

  my.p7 <- p7[[cancer.type]]
  stopifnot(ncol(my.p7) > 0)
  my.ex <- ex[[cancer.type]]
  stopifnot(ncol(my.ex) > 0)

  one.ex <- my.ex[ , sample.index, drop = FALSE]
  sample.id <- colnames(one.ex)
  one.spect <- my.p7[ , sample.id, drop = FALSE]
  one.ex <- one.ex[ , 1]

  if (isTRUE(out.dir)) {
    out.dir <- paste(gsub("::", "_", sample.id, fixed = TRUE),
                     mutation.type, sep = "_")
    message("out.dir is now ", out.dir)
  }

  rr <-
    OneMAPAssignTest(spect                   = one.spect,
                     reference.exp           = one.ex,
                     cancer.type             = cancer.type ,
                     mutation.type           = mutation.type,
                     exposure.mutation.type  = exposure.mutation.type,
                     max.mc.cores            = max.mc.cores,
                     max.level               = max.level,
                     out.dir                 = out.dir,
                     p.thresh                = p.thresh,
                     m.opts                  = m.opts,
                     max.presence.proportion = max.presence.proportion,
                     sigs.prop               = sigs.prop)
  return(rr)
}


#' Run one test of \code{\link{MAPAssignActivity1}}.
#'
#' @param spect A single spectrum.
#'
#' @param reference.exp Compare the inferred exposures to this.
#'
#' @param cancer.type Character string from a fixed set indicating
#'   different cancer types, used to look up the set of signatures
#'   known in that cancer type and the proportion of cancers of that
#'   type that have the signature. TODO: provide information on
#'   how to find the allowed cancer types.
#'
#' @param mutation.type One of "SBS96", "SBS192", "ID", "DBS78".
#'
#' @param exposure.mutation.type One of "SBS96", "ID", "DBS78".
#'
#' @param max.level The maximum number of signatures to try removing.
#'
#' @param out.dir If non-NULL create this directory if necessary and put
#'   results there.
#'
#' @param p.thresh If
#'  the p value for a better reconstruction with than without a set of signatures
#'  is > than \code{p.thresh}, then we can use exposures without this set.
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'
#' @param max.mc.cores
#'   The maximum number of cores to use.
#'   On Microsoft Windows machines it is silently changed to 1.
#'
#' @param max.subsets The maximum number of subsets that can be
#'   tested for removal from the set of signatures.
#'
#' @param sigs.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being signature identifiers. Can be the
#'    return value from \code{\link{ExposureProportions}}.
#'
#' @param sigs Matrix of signatures.
#'
#' @keywords internal
#'

OneMAPAssignTest <- function(spect,
                             reference.exp,
                             cancer.type,
                             mutation.type,
                             exposure.mutation.type,
                             max.subsets  = 1000,
                             max.level    = 5,
                             max.mc.cores = 100,
                             m.opts       = DefaultManyOpts(),
                             out.dir      = NULL,
                             p.thresh,
                             sigs.prop    = NULL,
                             sigs         = NULL) {
  if (!is.null(out.dir)) {
    if (!dir.exists(out.dir)) {
      created <- dir.create(out.dir, recursive = TRUE)
      stopifnot(created)
    }
    logfile <- file.path(out.dir, "log.txt")
    cat("date,\"", date(), "\"\n", file = logfile, sep = "")
  }

  shortlog <- function(tag, ...) {
    if (!is.null(out.dir)) {
      cat(tag,  "\n", sep = "", file = logfile, append = TRUE)
      utils::capture.output(print(...), file = logfile, append = TRUE)
      cat("\n", file = logfile, append = TRUE)
    }
  }

  if (is.null(sigs)) {
    sigs <- PCAWG7::signature$genome[[mutation.type]]
  }

  if (is.null(sigs.prop)) {
    sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                     cancer.type   = cancer.type,
                                     all.sigs = sigs)
  }

  mapout.time <- system.time(
    MAPout <-
      MAPAssignActivity1(
        spect                   = spect,
        sigs                    = sigs,
        sigs.presence.prop      = sigs.prop,
        max.level               = max.level, # length(sigs.prop) - 1,
        p.thresh                = p.thresh,
        m.opts                  = m.opts,
        max.mc.cores            = max.mc.cores, # mc.cores.per.sample = 100)
        max.subsets             = max.subsets)
  )
  shortlog("Time", mapout.time)

  if (!is.null(out.dir)) {
    save(MAPout, file = file.path(out.dir, "saved.MAPout.Rdata"))
  }

  if (!MAPout$success) {
    message("MAPAssignActivity1 did not succeed")
    if (!is.null(out.dir)) {
      cat("No result\n", MAPout$messages, "\n",
          file = file.path(out.dir, "no.results.txt"))
    }
    return(NULL)
  }

  if (is.matrix(reference.exp)) {
    if (ncol(reference.exp != 1)) {
      stop("Need 1-column matrix or vector for reference.exp")
    }
    reference.exp <- reference.exp[ , 1, drop = TRUE] # Make it a vector
  }
  ref.nonzero <-reference.exp[reference.exp > 0]
  tmp.names <- names(ref.nonzero)
  stopifnot(!is.null(tmp.names))
  ref.exp <- tibble::tibble(sig.id = tmp.names, ref.nonzero)
  rm(tmp.names)

  QP.exp <- OptimizeExposureQP(spect,
                               sigs[ , MAPout$MAP$sig.id,
                                     drop = FALSE])
  QP.best.MAP.exp <-
    tibble::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)

  qp.sparse <- OptimizeExposureQP(spect,
                                  sigs[ , MAPout$best.sparse$sig.id,
                                        drop = FALSE])
  QP.sparse.MAP.exp <-
    tibble::tibble(sig.id = names(qp.sparse), QP.sparse.MAP.exp = qp.sparse)

  MAP <- MAPout$MAP
  colnames(MAP)[2] <- "MAP.count"
  best.sparse <- MAPout$best.sparse
  colnames(best.sparse)[2] <- "best.sparse.count"

  comp <-
    dplyr::full_join(
      dplyr::full_join(
        dplyr::full_join(ref.exp, MAP),
        dplyr::full_join(best.sparse, QP.best.MAP.exp)),
      QP.sparse.MAP.exp)
  print(comp)

  if (!is.null(out.dir)) {
    out.path <- file.path(out.dir, "comparisions.csv")
    cat("\nComparisons of exposure attributions\n", file = out.path)
    suppressWarnings(
      utils::write.table(comp, file = out.path, append = TRUE, sep = ",",
                         row.names = FALSE))
  }

  # PCAWG attributions
  r.p <-
    ReconstructSpectrum(sigs, exp = ref.nonzero, use.sig.names = TRUE)

  # Best MAP
  r.b <- MAPout$MAP.recon
  # ReconstructSpectrum(sigs, exp = MAPout$MAP$count, use.sig.names = TRUE)

  # Best MAP re-optimized to a different objective function
  r.qp <-
    ReconstructSpectrum(sigs, exp = QP.exp, use.sig.names = TRUE)

  # MAP most sparse
  r.sparse.best <- MAPout$sparse.MAP.recon
  #  ReconstructSpectrum(
  #    sigs, exp = MAPout$best.sparse$count, use.sig.names = TRUE)

  # MAP most sparse re-optimized to different objective function
    r.qp.sparse <-
      ReconstructSpectrum(sigs, exp = qp.sparse, use.sig.names = TRUE)

  class(spect) <- "matrix"
  sol.matrix <- cbind(spect, r.p, r.b, r.qp, r.sparse.best, r.qp.sparse)

  colnames(sol.matrix) <-
    c("spect", "PCAWG7", "MAP",  "MAP+QP", "sparse", "sparse+QP")

  e.dist <- philentropy::distance(t(sol.matrix), use.row.names = TRUE)
  print(e.dist)
  if (!is.null(out.dir)) {
    cat("\nEuclidean distances\n", file = out.path, append = TRUE)
    suppressWarnings(
      utils::write.table(e.dist, file = out.path, append = TRUE, sep = ",",
                         col.names = NA, row.names = TRUE))
  }

  cos.sim <- philentropy::distance(t(sol.matrix),
                                   use.row.names = TRUE,
                                   method = "cosine")
  print(cos.sim)
  if (!is.null(out.dir)) {
    cat("\nCosine similarities\n", file = out.path, append = TRUE)
    suppressWarnings(
      utils::write.table(cos.sim, file = out.path, append = TRUE, sep = ",",
                         col.names = NA, row.names = TRUE))
    cat("Log likelihood of spectrum given reconstruction\n")
  }

  for (cname in colnames(sol.matrix)[2:ncol(sol.matrix)]) {
    logLH <-
      LLHSpectrumNegBinom(as.vector(spect),
                          sol.matrix[ , cname, drop = TRUE],
                          nbinom.size = m.opts$nbinom.size)
    cat("LogLH of spect from", cname, logLH, "\n")
    if (!is.null(out.dir)) {
      cat("spect", cname, logLH, "\n", sep = ",", file = out.path, append = TRUE)
    }
  }

  print(MAPout$MAP.row)
  print(MAPout$MAP.row$loglh.of.exp)

  colnames(sol.matrix) <- paste(colnames(sol.matrix), round(cos.sim[1, ], digits = 4))
  colnames(sol.matrix)[1] <- colnames(spect)
  if (!is.null(out.dir)) {
    tmp.catalog <- ICAMS::as.catalog(round(sol.matrix))
    ICAMS::PlotCatalogToPdf(
      tmp.catalog,
      file = file.path(out.dir, "reconstructions.pdf"))
    ICAMS::WriteCatalog(tmp.catalog,
                        file = file.path(out.dir, "reconstructions.csv" ))

  }

  print(MAPout$MAP.distances)
  print(MAPout$sparse.MAP.distances)

  return(MAPout)

}

