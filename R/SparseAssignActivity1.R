#' Component of \code{\link{SparseAssignActivity}} for one spectrum.
#' @keywords internal
#'
#' @param spect A single spectrum.
#'
#' @param sigs A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog.
#'
#' @param max.level The maximum number of signatures to try removing.
#'
#' @param p.thresh The p value threshold for deciding if a set of signatures is necessary.
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'
#' @param max.mc.cores
#'   The maximum number of cores to use.
#'   On Microsoft Windows machines it is silently changed to 1.)


SparseAssignActivity1 <- function(spect,
                                  sigs,
                                  max.level    = 5,
                                  p.thresh     = 0.05,
                                  m.opts,
                                  max.mc.cores = min(20, 2^max.level)) {

  my.msg <- function(trace.level, ...)
    if (m.opts$trace >= trace.level) message("SparseAssignActivity1: ", ...)

  max.sig.index <- ncol(sigs)
  my.msg(10, "SparseAssignActivity1: number of signatures = ", max.sig.index)
  mode(spect) <-  'numeric'
  start <- OptimizeExposure(spect,
                            sigs,
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
  non.0.exp.index <- which(start.exp > 0.5) # indices of signatures with non-zero
                                            # exposures； TODO, possibly remove
  my.msg(0, "Starting with ",
         paste(names(start.exp)[non.0.exp.index], collapse = ","),
         "\nmax.level = ", max.level,
         "\nlog likelihood using all signatures = ", lh.w.all)

  if (length(non.0.exp.index) < 2) {
    my.msg(0, "returning, only ", length(non.0.exp.index), " non-0 exposures")
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)

  mc.cores <- Adj.mc.cores(max.mc.cores) # Set to 1 if OS is MS Windows

  # c.r.s is "cannot remove subsets", i.e. subset of the signature indices that
  # cannot be removed
  c.r.s <- sets::set()
  best.exp <- list() # Index by level, the expression of the assignments (exposures) at a given level
  best.sig.indices <- list() # Index by level, the signature indices for the exposures in best.exp
  # TODO, Should best actually be maximum likelihood or maximum a posteriori?

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

    # Get the maximum likelihood exposure for try.sigs
    try.exp <- OptimizeExposure(spect,
                                try.sigs,
                                m.opts = m.opts)

    # TODO -- deal with case when try$loglh is Inf

    statistic <- 2 * (lh.w.all - try.exp$loglh)
    chisq.p <- stats::pchisq(statistic, df, lower.tail = FALSE)
    return(list(p                 = chisq.p,
                exp               = try.exp[["exposure"]],
                sig.indices       = try.sigs.indices,
                removed.sig.names = colnames(sigs)[subset.to.remove.v],
                loglh.of.exp      = try.exp[["loglh"]]
                ))
  }

  for (df in 1:max.level) {
    my.msg(0, "\ndf = ", df)
    subsets <- as.list(sets::set_combn(non.0.exp.index, df))
    discard <- lapply(subsets, is.superset.of.any, background = c.r.s)
    subsets2 <- subsets[!(unlist(discard))]
    if (length(subsets2) == 0) break;

    my.msg(0, "Number of subsets to remove = ", length(subsets2))
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

    p.to.remove <- unlist(lapply(check.to.remove, `[`, "p"))
    names(p.to.remove) <-
      unlist(lapply(check.to.remove,
                    function(x) paste(x$removed.sig.names, collapse = ",")))
    if (m.opts$trace > 10) {
      message("SparseAssignActivity1: p.to.remove =")
      for (ii in 1:length(p.to.remove)) {
        message(ii, " ", names(p.to.remove)[ii], " ", p.to.remove[ii])
      }
    }
    if (all(p.to.remove < p.thresh)) {
      my.msg(10, "Cannot remove any subsets at level ", df)
      break;
    }
    cannot.remove <- subsets2[p.to.remove < p.thresh]
    c.r.s <- sets::set_union(c.r.s, sets::as.set(cannot.remove))
    xx <- which(p.to.remove == max(p.to.remove))
    my.msg(10, "xx = ", paste(xx, collapse = ","),
           "\nmax(p.to.remove = ", max(p.to.remove), ")")
    if (length(xx) > 1) {
      possible.sigs <-
        lapply(X = subsets[xx],
               function(subset) paste(colnames(sigs)[unlist(subset)], collapse = ","))
      xx <- min(xx)
      sig.name.to.remove <- colnames(sigs)[unlist(subsets[[xx]])]

      message("Temporary warning, > 1 signature can be removed; ",
              "selecting the first one arbitrarily (index = ", xx, ",
              name = ", paste(sig.name.to.remove, collapse = ","), ")\n",
              "Possibilities were\ ", paste(possible.sigs, collapse = "\n"))
    }
    best.exp[df] <- list(check.to.remove[[xx]]$exp)
    best.sig.indices[df] <- list(check.to.remove[[xx]]$sig.indices)
  }

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
  stopifnot(abs(sum(out.exp) - sum(spect)) < 1)
  return(out.exp)
}
