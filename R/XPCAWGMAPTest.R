#' Run \code{\link{MAPAssignActivity1}} on one sample from the PCAWG platinum data set with and with pre-filtering by bootstrapped quadratic programming.
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
#' @return A list with two elements, each the result for one call to \code{\link{OneMAPAssignTest}}.
#'
#' @export
#'

XPCAWGMAPTest <- function(cancer.type,
                         sample.index,
                         mutation.type,
                         max.level = 5,
                         max.mc.cores,
                         out.dir = NULL,
                         p.thresh = 0.01,
                         m.opts = DefaultManyOpts(),
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

  if (sample.index > ncol(my.ex)) {
    message("Sample index (", sample.index, ") > ncol(my.ex) (", ncol(my.ex), ")")
    message("Returning NULL")
    return(NULL)
  }
  one.ex <- my.ex[ , sample.index, drop = FALSE]
  sample.id <- colnames(one.ex)
  one.spect <- my.p7[ , sample.id, drop = FALSE]
  one.ex <- one.ex[ , 1]

  if (isTRUE(out.dir)) {
    out.dir <- paste(gsub("::", "_", sample.id, fixed = TRUE),
                     mutation.type, sep = "_")
    message("out.dir is now ", out.dir)
  }

  sigs <- PCAWG7::signature$genome[[mutation.type]]

  if (is.null(sigs.prop)) {
    sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                     cancer.type   = cancer.type,
                                     all.sigs = sigs)
  }

  if (TRUE) {
    # qp.assign <- OptimizeExposureQP(spect, sigs) Not used a this point
    qp.assign <-
      ICAMS.shiny::GetExposureWithConfidence(
        catalog = one.spect,
        sig.universe = sigs[ , names(sigs.prop), drop = FALSE],
        conf.int = 0.75)
    some.sigs.prop <- sigs.prop[rownames(qp.assign)]
  } else {
    some.sigs.prop <- sigs.prop[setdiff(names(sigs.prop), PossibleArtifacts())]
  }

  one.assign <- function(my.sigs.prop, subdir) {
    rr <-
      OneMAPAssignTest(spect                   = one.spect,
                       reference.exp           = one.ex,
                       cancer.type             = cancer.type ,
                       mutation.type           = mutation.type,
                       exposure.mutation.type  = exposure.mutation.type,
                       max.mc.cores            = max.mc.cores,
                       max.level               = max.level,
                       out.dir                 = file.path(out.dir, subdir),
                       p.thresh                = p.thresh,
                       m.opts                  = m.opts,
                       max.presence.proportion = max.presence.proportion,
                       sigs.prop               = my.sigs.prop)
    return(rr)
  }

  message("Testing all sigs ", paste(names(sigs.prop), collapse = ", "))
  rr1 <- one.assign(sigs.prop, "all.sigs")

  if (length(some.sigs.prop) == length(sigs.prop)) {
    message("No signatures removed")
    rr2 <- NULL
  } else {
    message("Testing ", paste(names(some.sigs.prop), collapse = ", "))
    rr2 <- one.assign(some.sigs.prop, "some.sigs")
    if (!isTRUE(all.equal(rr1$MAP, rr2$MAP)) && FALSE) {
       for (sig.name in setdiff(names(sigs.prop), names(some.sigs.prop))) {
            # Later try to salvage by a forward search?

       }

    }

  }

  if (!is.null(out.dir) && !is.null(rr2)) {
    if (!isTRUE(all.equal(rr1$MAP, rr2$MAP))) {
      cat("DIFFERNCE.MAP\n",
          file = file.path(out.dir, "DIFFERENCE.MAP"))
    } else {
      cat("OK.MAP\n",
          file = file.path(out.dir, "OK.MAP"))
    }
    if (!isTRUE(all.equal(rr1$best.sparse, rr2$best.sparse))) {
      cat("DIFFERNCE.best.sparse\n",
          file = file.path(out.dir, "DIFFERENCE.best.sparse"))
    } else {
      cat("OK.best.sparse\n",
          file = file.path(out.dir, "OK.best.sparse"))
    }
  }

  return(list(all = rr1, some = rr2))
}
