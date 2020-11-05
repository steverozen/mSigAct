#' Component of \code{\link{SparseAssignActivity}} for one spectrum.
#' @keywords internal
#'
#' @param spect A single spectrum.
#'
#' @param sigs A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog.
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
#'   If \code{NULL} defaults to \code{2^max.level} -- except on
#'    MS Windows machines, where it defaults to 1.)


SparseAssignActivity1 <- function(spect,
                                  sigs,
                                  max.level    = 5,
                                  p.thresh     = 0.05,
                                  eval_f,
                                  m.opts,
                                  max.mc.cores = NULL) {
  mode(spect) <-  'numeric'
  start <- OptimizeExposure(spect,
                            sigs,
                            eval_f  = eval_f,
                            m.opts  = m.opts)
  lh.w.all <- start$loglh  # The likelihood with all signatures
  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp > 0.5) # indices of signatures with non-zero
                                            # expsoures； TODO, possibly remove
  m.opts$trace <- 100
  if (m.opts$trace > 0) {
    message('Starting with ',
            paste(names(start.exp)[non.0.exp.index], collapse = ","),
            '\n')
    message('max.level =', max.level, '\n')
    message("full -log likelihood is ", lh.w.all)
  }
  if (length(non.0.exp.index) < 2) {
    if (m.opts$trace > 0)
      message("returning, only ", length(non.0.exp.index),
              " non-0 exposures")
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)

  if (is.null(max.mc.cores)) { max.mc.cores <- 2^max.level }

  mc.cores <- Adj.mc.cores(max.mc.cores) # Set to 1 if OS is MS Windows

  c.r.s <- sets::set() # subsets of the signature indices that cannot be removed

  best.exp <- list()
  best.sig.indices <- list()

  info.of.removed.subset <- function(subset) {
    # subset is of class set (package sets)
    subset.to.remove.v <- as.numeric(subset)
    tmp.set <- setdiff(non.0.exp.index, subset.to.remove.v)
    try.sigs <- sigs[ , tmp.set, drop = FALSE]
    if (m.opts$trace > 0) {
      message("\nTesting removal of signature index(s) ",
              paste(subset.to.remove.v, collapse = ", "),
              "; name(s) = ",
              paste(colnames(sigs)[subset.to.remove.v], collapse = ", "))
    }

    # Get the max lh for try.sig
    try <- OptimizeExposure(spect,
                          try.sigs,
                          eval_f = eval_f,
                          m.opts = m.opts)

    # TODO -- deal with case when try$loglh is Inf

    statistic <- 2 * (lh.w.all - try$loglh)
    chisq.p <- stats::pchisq(statistic, df, lower.tail = FALSE)
    return(list(p = chisq.p, exp = try$exposure, sig.indices = tmp.set))
  }

  for (df in 1:max.level) {
    if (m.opts$trace > 0) message("\ndf = ", df)
    subsets <- as.list(sets::set_combn(non.0.exp.index, df))
    discard <- lapply(subsets, is.superset.of.any, background = c.r.s)
    subsets2 <- subsets[!(unlist(discard))]
    if (length(subsets2) == 0) break;

    check.to.remove <-
      parallel::mclapply(X = subsets2,
                         FUN = info.of.removed.subset,
                         mc.cores = mc.cores)

    p.to.remove <- unlist(lapply(check.to.remove, function(x) x$p))

    if (all(p.to.remove < p.thresh)) break;
    cannot.remove <- subsets2[p.to.remove < p.thresh]
    c.r.s <- sets::set_union(c.r.s, sets::as.set(cannot.remove))
    xx <- which(p.to.remove == max(p.to.remove))
    if (length(xx) > 1) {
      xx <- min(xx)
      sig.name.to.remove <- colnames(sigs)[xx]
      warning("Temporary warning, > 1 signature can be removed; ",
              "selecting the first one arbitrarily (index = ", xx, ",
              name = ", sig.name.to.remove)
    }
    best.exp[df] <- list(check.to.remove[[xx]]$exp)
    best.sig.indices[df] <- list(check.to.remove[[xx]]$sig.indices)
  }

  # Need to return the exposures in the context of the
  # orginal signatures matrix
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


#' Is the set 'probe' a superset of any set in 'background'?
#'
#' @keywords internal
is.superset.of.any <- function(probe, background) {
  for (b in background) {
    if (sets::set_is_proper_subset(b, probe)) return(TRUE)
  }
  return(FALSE)
}
