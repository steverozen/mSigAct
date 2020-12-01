#' Run \code{\link{MAPAssignActivity1}} on one sample from the PCAWG platinum data set with global opts maxeval 10000 and 1000 and compare the results
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
#' @param eval_f See \code{\link[nloptr]{nloptr}}.
#'
#' @param eval_g_ineq See \code{\link[nloptr]{nloptr}}.
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

YPCAWGMAPTest <- function(cancer.type,
                         sample.index,
                         mutation.type,
                         max.level = 5,
                         max.mc.cores,
                         out.dir = NULL,
                         p.thresh = 0.01,
                         m.opts = DefaultManyOpts(),
                         eval_f = ObjFnBinomMaxLHRound,
                         eval_g_ineq = NULL,
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

  one.assign <- function(my.sigs.prop, gme) {
    my.m.opts <- m.opts
    my.m.opts$global.opts$maxeval <- gme
    message("my.m.opts$global.opts$maxeval = ", gme)
    my.out.dir <- file.path(out.dir, paste0("d", gme))
    rr <-
      OneMAPAssignTest(spect                   = one.spect,
                       reference.exp           = one.ex,
                       cancer.type             = cancer.type ,
                       mutation.type           = mutation.type,
                       exposure.mutation.type  = exposure.mutation.type,
                       max.mc.cores            = max.mc.cores,
                       max.level               = max.level,
                       out.dir                 = my.out.dir,
                       p.thresh                = p.thresh,
                       m.opts                  = my.m.opts,
                       eval_f                  = eval_f,
                       eval_g_ineq             = eval_g_ineq,
                       max.presence.proportion = max.presence.proportion,
                       sigs.prop               = my.sigs.prop)
    print(rr$time.for.MAP.assign)
    if (!is.null(out.dir)) {
      cat(rr$time.for.MAP.assign, "\n", file = file.path(my.out.dir, "timing.txt"))
    }
    return(rr)
  }

  rr1 <- one.assign(sigs.prop, 10000)

  rr2 <- one.assign(sigs.prop, 1000)

  if (!is.null(out.dir) && !is.null(rr2)) {
    if (!isTRUE(all.equal(rr1$MAP, rr2$MAP, tol = 3))) {
      cat("DIFFERNCE.MAP\n",
          file = file.path(out.dir, "DIFFERENCE.MAP"))
    } else {
      cat("OK.MAP\n",
          file = file.path(out.dir, "OK.MAP"))
    }
    if (!isTRUE(all.equal(rr1$best.sparse, rr2$best.sparse, tol = 3))) {
      cat("DIFFERNCE.best.sparse\n",
          file = file.path(out.dir, "DIFFERENCE.best.sparse"))
    } else {
      cat("OK.best.sparse\n",
          file = file.path(out.dir, "OK.best.sparse"))
    }
  }

  return(list(all = rr1, some = rr2))
}
