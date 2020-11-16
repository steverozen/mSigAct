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
#' @param out.dir If non-NULL creat this directory if necessary and put
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
#' @param max.presence.proportion xxxxx
#'
#' @export
#'

PCAWGMAPTest <- function(cancer.type,
                         sample.index,
                         mutation.type,
                         max.level = 5,
                         max.mc.cores,
                         out.dir = NULL,
                         p.thresh = 0.01,
                         m.opts = DefaultManyOpts(),
                         max.presence.proportion = 0.99) {

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
                     max.presence.proportion = max.presence.proportion)
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
#' @param out.dir If non-NULL creat this directory if necessary and put
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
#' @param max.subsets The maxium number of subsets that can be
#'   tested for removal from the set of signatures.
#'
#' @param max.presence.proportion The maxium value of the proportion
#'   of tumors that must have a given signature.
#'
#' @export
#'
OneMAPAssignTest <- function(spect,
                             reference.exp,
                             cancer.type,
                             mutation.type,
                             exposure.mutation.type,
                             max.subsets = 1000,
                             max.level   = 5,
                             max.mc.cores = 100,
                             m.opts = DefaultManyOpts(),
                             out.dir = NULL,
                             p.thresh,
                             max.presence.proportion) {

  if (!is.null(out.dir)) {
    if (!dir.exists(out.dir)) {
      created <- dir.create(out.dir, recursive = TRUE)
      stopifnot(created)
    }
  }

  sigs.prop <- PCAWG7::exposure.stats$PCAWG[[exposure.mutation.type]][[cancer.type]]
  if (is.null(sigs.prop)) {
    stop("Cannot find exposure information for ",
         exposure.mutation.type, " for ", cancer.type)
  }
  sig.names <- rownames(sigs.prop)
  sigs.prop <- unlist(sigs.prop[ , 2])
  names(sigs.prop) <- sig.names

  sigs.no.info <-
    setdiff(sig.names,
            colnames(PCAWG7::signature$genome[[mutation.type]]))

  for (zz in sigs.no.info) {
    message("No signature ", zz, " for ", mutation.type)
    message("Dropping from signature universe")
    sigs.prop <- sigs.prop[-(which(names(sigs.prop) == zz))]
  }

  sigs <- PCAWG7::signature$genome[[mutation.type]][ , names(sigs.prop), drop = FALSE]

  qp.assign <- OptimizeExposureQP(spect, sigs)

  MAPout <-
    mSigAct::MAPAssignActivity1(
      spect = spect,
      sigs = sigs,
      sigs.presence.prop = sigs.prop,
      max.level = max.level, # length(sigs.prop) - 1,
      p.thresh = p.thresh,
      eval_f = ObjFnBinomMaxLHNoRoundOK,
      m.opts = m.opts,
      max.mc.cores = max.mc.cores, # mc.cores.per.sample = 100)
      max.subsets = max.subsets,
      max.presence.proportion = max.presence.proportion)

  if (is.null(MAPout)) {
    message("No result from MAPAssignActivity1")
    if (!is.null(out.dir)) {
      cat("No result\n", file = file.path(out.dir, "no.results.txt"))
    }
    return(NULL)
  }

  xx <- ListOfList2Tibble(MAPout)

  if (!is.null(out.dir)) {
    saved.MAPout.list    <- MAPout
    saved.MAPout.tribble <- xx
    save(saved.MAPout.list, saved.MAPout.tribble,
         file = file.path(out.dir, "all.saved.results.Rdata"))
  }

  ref.nonzero <-reference.exp[reference.exp > 0]
  ref.exp <- tibble::tibble(sig.id = names(ref.nonzero), ref.nonzero)

  # select.best and todo compare to PCAWG exposure
  best <- dplyr::arrange(xx, rlang::.data$MAP)[nrow(xx),  ]
  names.best <- names(best[["exp"]])
  best.exp <- best[["exp"]][[1]]
  if (is.null(names(best.exp))) {
    names(best.exp) <- names.best
  }
  MAP.best.exp <- tibble::tibble(sig.id = names(best.exp), best.exp )

  sparse.best <-
    dplyr::arrange(xx, rlang::.data$df, rlang::.data$MAP)[nrow(xx), ]
  names.sparse.best <- names(sparse.best[["exp"]])
  most.sparse.exp <- sparse.best[["exp"]][[1]]
  if (is.null(names(most.sparse.exp))) {
    names(most.sparse.exp) <- names.sparse.best # Necessary if only 1 signature
  }
  MAP.most.sparse <-
    tibble::tibble(sig.id = names(most.sparse.exp), most.sparse = most.sparse.exp)

  QP.exp <- OptimizeExposureQP(spect, sigs[ , names(best.exp), drop = FALSE])
  QP.best.MAP.exp <-
    tibble::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)

  qp.sparse <- OptimizeExposureQP(spect, sigs[ , names(most.sparse.exp), drop = FALSE])
  QP.sparse.MAP.exp <-
    tibble::tibble(sig.id = names(qp.sparse), QP.sparse.MAP.exp = qp.sparse)

  comp <-
    dplyr::full_join(
      dplyr::full_join(
        dplyr::full_join(ref.exp, MAP.best.exp),
        dplyr::full_join(MAP.most.sparse, QP.best.MAP.exp)),
      QP.sparse.MAP.exp)
  print(comp)

  if (!is.null(out.dir)) {
    out.path <- file.path(out.dir, "comparisions.csv")
    cat("\nComparisons of exposure attributions\n", file = out.path)
    suppressWarnings(
      utils::write.table(comp, file = out.path, append = TRUE, sep = ","))
  }

  # PCAWG attributions
  r.p <-
    ReconstructSpectrum(sigs, exp = ref.nonzero, use.sig.names = TRUE)

  # Best MAP
  r.b <-
    ReconstructSpectrum(sigs, exp = best.exp, use.sig.names = TRUE)

  # Best MAP re-optimized to a different obective function
  r.qp <-
    ReconstructSpectrum(sigs, exp = QP.exp, use.sig.names = TRUE)

  # MAP most sparse
  r.sparse.best <-
    ReconstructSpectrum(sigs, exp = most.sparse.exp, use.sig.names = TRUE)

  # MAP most sparse re-optimsed to different objective function
    r.qp.sparse <-
      ReconstructSpectrum(sigs, exp = qp.sparse, use.sig.names = TRUE)

  sol.matrix <- cbind(spect, r.p, r.b, r.qp, r.sparse.best, r.qp.sparse)

  colnames(sol.matrix) <-
    c("spect", "PCAWG7", "MAP",  "MAP+QP", "sparse", "sparse+QP")

  e.dist <- philentropy::distance(t(sol.matrix), use.row.names = TRUE)
  print(e.dist)
  if (!is.null(out.dir)) {
    cat("\nEuclidean distances\n", file = out.path, append = TRUE)
    suppressWarnings(
      utils::write.table(e.dist, file = out.path, append = TRUE, sep = ","))
  }

  cos.sim <- philentropy::distance(t(sol.matrix),
                                   use.row.names = TRUE,
                                   method = "cosine")
  print(cos.sim)
  if (!is.null(out.dir)) {
    cat("\nCosine similarities\n", file = out.path, append = TRUE)
    suppressWarnings(
      utils::write.table(cos.sim, file = out.path, append = TRUE, sep = ","))
  }

  colnames(sol.matrix) <- paste(colnames(sol.matrix), round(cos.sim[1, ], digits = 4))
  colnames(sol.matrix)[1] <- colnames(spect)
  if (!is.null(out.dir)) {
    ICAMS::PlotCatalogToPdf(
      ICAMS::as.catalog(round(sol.matrix)),
      file = file.path(out.dir, "reconstructions.pdf"))
  }

}

