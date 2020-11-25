#' Optimize assignment of a fixed set of signature activities for a \strong{single} tumor.
#'
#' @inheritParams OptimizeExposure
#'
#' @keywords internal
Nloptr1Tumor <- function(spectrum,
                         sigs,
                         m.opts = NULL,
                         eval_f,
                         ... ) {
  if (!"matrix" %in% class(sigs)) {
    if (!"numeric" %in% class(sigs)) {
      stop("Unexpected class argument: ", class(sigs))
    }
    # Otherwise sigs is a numeric vector, not a matrix. We assume the caller
    # intended it as a single column matrix.
    sigs <- matrix(sigs, ncol = 1)
  }

  if (nrow(sigs) != length(spectrum)) {
    # If we get here there is an error. We save spectrum and sigs, which seems
    # useful for debugging in a call to mclapply (multi-core lapply).
    save(spectrum, sigs, file = 'spectrum.sigs.debug.Rdata')
    stop("nrow(sigs) != length(spectrum), look in file spectrum.sigs.debug.Rdata")
  }

  stopifnot(mode(spectrum) == "numeric")
  foo <- spectrum
  storage.mode(spectrum) <- "double"
  stopifnot(foo == spectrum)

  if (is.null(m.opts)) m.opts <- DefaultManyOpts()

  use.old <- TRUE
  if (ncol(sigs) == 1 || use.old) {
    # x0 is uniform distribution of mutations across signatures
    # (Each signature gets the same number of mutations)
    x0 <- rep(sum(spectrum) / ncol(sigs), ncol(sigs))
  } else {
    stop("This branch does not work")
    # We get x0 that minimizes the Frobenius norm
    # x0 <- sum(spectrum) * as.vector(ICAMS.shiny:::findSigExposures(M = spectrum, P = sigs)$exposures)
    message("QP estimate is ", paste(x0, collapse = " "))
    m.opts$global.opts$maxeval <- 1
    m.opts$local.opts$maxeval <- 30000  }

  if (!is.null(m.opts$global.opts)) { # WARNING, ADDITIONAL CODE WOULD NEED TO BE CHANGED DISABLE THIS BRANCH

    global.res <- nloptr::nloptr(
      x0       = x0,
      eval_f   = eval_f,
      lb       = rep(0, ncol(sigs)),
      ub       = rep(sum(spectrum), ncol(sigs)),
      opts     = m.opts$global.opts,
      spectrum = spectrum,
      sigs     = sigs,
      ...)
    if (m.opts$trace > 0) {
      message("global.res$objective = ", global.res$objective)
    }
    if (global.res$iterations == m.opts$global.opts[["maxeval"]]) {
      # warning("reached maxeval on global optimization: ", global.res$iterations)
    }
  }

  if (use.old) {
    my.x0 <- global.res$solution
  } else {
    my.x0 <- x0
  }
  local.res <- nloptr::nloptr(
    x0=my.x0,
    # x0       = global.res$solution,
    eval_f   = eval_f,
    lb       = rep(0, ncol(sigs)),
    ub       = rep(sum(spectrum) + 1e-2, ncol(sigs)),
    opts     = m.opts$local.opts,
    spectrum = spectrum,
    sigs     = sigs,
    ...)
  if (m.opts$trace > 0)
    message("local.res$objective = ", local.res$objective)
  message("local.res$iterations = ", local.res$iterations)

  maxeval.warning <- NULL
  if (local.res$iterations == m.opts$local.opts[["maxeval"]]) {
    if (m.opts$trace > 0)
      message("reached maxeval on local optimization: ", local.res$iterations)
    maxeval.warning <- list(spectrum = spectrum, sigs = sigs)
  }

  names(local.res$solution) <- colnames(sigs)
  return(list(objective      =  local.res$objective,
              solution       = local.res$solution,
              global.res     = global.res,
              local.res      = local.res,
              maxeval.warnng = maxeval.warning))
}

