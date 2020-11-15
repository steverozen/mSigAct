#' Component of \code{\link{SparseAssignActivity}} for one spectrum.
#'
#' @export
#'
#' @param spect A single spectrum.
#'
#' @param sigs A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being the same as \code{colnames(sigs)}.
#'
#' @param max.level The maximum number of signatures to try removing.
#'
#' @param eval_f See \code{\link[nloptr]{nloptr}}.
#'
#' @param p.thresh The p value threshold for deciding if a set of signatures is necessary.
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'
#' @param max.mc.cores
#'   The maximum number of cores to use.
#'   On Microsoft Windows machines it is silently changed to 1.)


MAPAssignActivity1 <- function(spect,
                               sigs,
                               sigs.presence.prop,
                               max.level    = 5,
                               p.thresh     = 0.05,
                               eval_f,
                               m.opts,
                               max.mc.cores = min(20, 2^max.level),
                               max.subsets = 1000) {

  my.msg <- function(trace.level, ...)
    if (m.opts$trace >= trace.level) message("MAPAssignActivity1: ", ...)

  max.sig.index <- ncol(sigs)
  my.msg(10, "number of signatures = ", max.sig.index)
  mode(spect) <-  'numeric'
  stopifnot(isTRUE(all.equal(colnames(sigs), names(sigs.presence.prop))))
  start <- OptimizeExposure(spect,
                            sigs,
                            eval_f  = eval_f,
                            m.opts  = m.opts)
  lh.w.all <- start$loglh  # The likelihood with all signatures
  if (lh.w.all == -Inf) {
    problem <- "Cannot generate spect from sigs"
    message(problem)
    rr <- rep(NaN, ncol(sigs))
    attr(rr, "log.likelihood") <- lh.w.all
    attr(rr, "problem") <- problem
    return(rr)
  }
  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp >= 0.5) # indices of signatures with non-zero
                                             # exposures

  if (length(non.0.exp.index) < ncol(sigs)) {
    removed.sig.names <- colnames(sigs)[which(start.exp < 0.5)]
    my.msg(0, ncol(sigs) - length(non.0.exp.index),
           " signatures removed at beginning")
    my.msg(0, "removed: ",
           paste(removed.sig.names, collapse = ", "))
  }

  df0.sig.names <- colnames(sigs)[non.0.exp.index]
  my.msg(0, "Starting with ",
         paste(df0.sig.names, collapse = ","),
         "\nmax.level = ", max.level,
         "\nlog likelihood using all signatures = ", lh.w.all)

  df0.prob.of.model <- P.of.M(df0.sig.names, sigs.presence.prop)

  df0 <- list(sig.names            = paste(df0.sig.names, collapse = ","),
              p.for.removing.sigs  = NA,
              exp                  = start.exp[non.0.exp.index],
              sig.indices          = non.0.exp.index,
              removed.sig.names    = "",
              loglh.of.exp         = lh.w.all,
              prob.of.model        = df0.prob.of.model,
              MAP                  = lh.w.all + df0.prob.of.model,
              df                   = 0)


  out.list <- list(df0)

  if (length(non.0.exp.index) < 2) {
    my.msg(0, "returning, only ", length(non.0.exp.index), " non-0 exposures")
    return(out.list)
  }

  max.level <- min(max.level, length(non.0.exp.index) - 1)

  mc.cores <- Adj.mc.cores(max.mc.cores) # Set to 1 if OS is MS Windows

  # c.r.s is "cannot remove subsets", i.e. subset of the signature indices that
  # cannot be removed
  c.r.s <- sets::set()

  old.sparse <- FALSE

  if (old.sparse) {
    best.exp <- list() # Index by level, the expression of the assignments (exposures) at a given level
    best.sig.indices <- list() # Index by level, the signature indices for the exposures in best.exp
  }

  # Internal function ====================================================
  info.of.removed.subset <- function(subset) {
    # subset is of class set (package sets)
    subset.to.remove.v <- as.numeric(subset)
    if (any(subset.to.remove.v > max.sig.index)) {
      stop("Got illegal signature index in subset.to.remove = ",
           paste(subset.to.remove.v, collapse = " "))
    }
    try.sigs.indices <- setdiff(non.0.exp.index, subset.to.remove.v)
    try.sigs <- sigs[ , try.sigs.indices, drop = FALSE]

    my.msg(10, "Testing removal of signature index(s) ",
           paste(subset.to.remove.v, collapse = ", "),
           "; name(s) = ",
           paste(colnames(sigs)[subset.to.remove.v], collapse = ", "))

    try.sigs.names <- colnames(sigs)[try.sigs.indices]

    # Get the maximum likelihood exposure for try.sigs
    try.exp <- OptimizeExposure(spect,
                                try.sigs,
                                eval_f = eval_f,
                                m.opts = m.opts)

    if (is.infinite(try.exp$loglh)) {
      # It was not possible to generate the spectrum from the signatures, so
      # try.exp$loglh should be -Inf and the p value the test that the spectrum
      # can be better reconstructed with union(try.sigs, subset.to.remove.v)
      # than with try.sigs is 0.
       if (try.exp$loglh < 0) {
          chisq.p <- 0
       } else {
         stop("How did we here?")
       }
    } else {
      statistic <- 2 * (lh.w.all - try.exp$loglh)
      chisq.p <- stats::pchisq(statistic, df, lower.tail = FALSE)
    }
    prob.of.model <- P.of.M(try.sigs.names, sigs.presence.prop)

    return(list(
      sig.names            = paste(colnames(sigs)[try.sigs.indices],
                                   collapse = ","),
      p.for.removing.sigs  = chisq.p,
      exp                  = try.exp[["exposure"]],
      sig.indices          = try.sigs.indices,
      removed.sig.names    = paste(colnames(sigs)[subset.to.remove.v],
                                   collapse = ","),
      loglh.of.exp         = try.exp[["loglh"]],
      prob.of.model        = prob.of.model,
      MAP                  = try.exp[["loglh"]] + prob.of.model,
      df                   = df))
  } # Internal function info.of.removed.subset ===========================


  for (df in 1:max.level) {
    my.msg(0, "\ndf = ", df)
    subsets <- as.list(sets::set_combn(non.0.exp.index, df))
    discard <- lapply(subsets, is.superset.of.any, background = c.r.s)
    subsets2 <- subsets[!(unlist(discard))]
    if (length(subsets2) == 0) break;

    my.msg(0, "Number of subsets to remove = ", length(subsets2))
    if (length(subsets2) > max.subsets) {
      my.msg(-1, "Number of subsets (", length(subsets2), ") > max.subsets (", max.subsets, ")")
      my.msg(-1, "Returning NULL")
      return(NULL)
    }
    time.used <- system.time(
      check.to.remove <-
        parallel::mclapply(X = subsets2,
                           FUN = info.of.removed.subset,
                           mc.cores = max.mc.cores)
    )
    if (m.opts$trace > 3) {
      message("Time used to compute p values for df ", df, ":")
      print(time.used)
    }

    check.mclapply.result(check.to.remove, "SparseAssignActivity1")

    p.to.remove <- unlist(lapply(check.to.remove, `[`, "p.for.removing.sigs"))
    names(p.to.remove) <-
      unlist(lapply(check.to.remove,
                    function(x) paste(x$removed.sig.names, collapse = ",")))
    if (m.opts$trace > 10) {
      message("MAPAssignActivity1: p.to.remove = ")
      for (ii in 1:length(p.to.remove)) {
        message(ii, " ", names(p.to.remove)[ii], " ", p.to.remove[ii])
      }
    }
    if (all(p.to.remove < p.thresh)) {
      my.msg(10, "Cannot remove any subsets at level ", df)
      break;
    }
    ok.after.removal <- check.to.remove[p.to.remove >= p.thresh]
    out.list <-c(out.list, ok.after.removal)

    cannot.remove <- subsets2[p.to.remove < p.thresh]
    c.r.s <- sets::set_union(c.r.s, sets::as.set(cannot.remove))

    if (old.sparse) {
      # We keep track of the subset of signatures that are least important
      xx <- which(p.to.remove == max(p.to.remove))
      my.msg(10, "xx = ", paste(xx, collapse = ","),
             "\nmax(p.to.remove = ", max(p.to.remove), ")")

      if (length(xx) > 1) {
        possible.sigs <-
          lapply(X = subsets[xx],
                 function(subset) paste(colnames(sigs)[unlist(subset)], collapse = ","))
        xx <- min(xx)
        sig.name.to.remove <- colnames(sigs)[unlist(subsets[[xx]])]

        message("> 1 subset of signatures can be removed; ",
                "selecting the first one arbitrarily (index = ", xx,
                ",\nsigs are ", paste(sig.name.to.remove, collapse = ","), ")\n",
                "Possibilities were\ ", paste(possible.sigs, collapse = "\n"))
      }
      best.exp[df] <- list(check.to.remove[[xx]]$exp)
      best.sig.indices[df] <- list(check.to.remove[[xx]]$sig.indices)
    }

  }

  if (old.sparse) {
    # Need to return the exposures in the context of the
    # original signatures matrix
    out.exp <- numeric(ncol(sigs)) #all zeros
    names(out.exp) <- colnames(sigs)
    max.df <- length(best.sig.indices)
    if (max.df == 0) {
      out.exp <- start.exp
    } else {
      out.exp[unlist(best.sig.indices[[max.df]])] <- unlist(best.exp[[max.df]])
    }
    # stopifnot(abs(sum(out.exp) - sum(spect)) < 1)
    my.msg(0, "sum(out.exp) - sum(spect) = ", sum(out.exp) - sum(spect))
  }

  return(out.list)
}

#' Calculate \eqn{P(M)} -- the probabily of a model of which signatures are present in a sample.
#'
#' @param model Names of sigs present in a trial exposure. Do not use indices.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being a superset of \code{model}.
#'
#' @keywords internal
#'
P.of.M <- function(model, sigs.presence.prop) {
  stopifnot(length(setdiff(model, names(sigs.presence.prop))) == 0)
  present <- sigs.presence.prop[model]
  not.present.names <- setdiff(names(sigs.presence.prop), model)
  not.present <- sigs.presence.prop[not.present.names]
  not.present <- 1 - not.present
  rr <- sum(log(c(present, not.present)))
  # This could be -Inf if one of the signatures is always present in the prior data but absent; probably want to fudge this slightly
  return(rr)
}
