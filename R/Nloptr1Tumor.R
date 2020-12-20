#' Optimize assignment of a fixed set of signature activities for a \strong{single} tumor.
#'
#' @inheritParams OptimizeExposure
#'
#' @keywords internal
Nloptr1Tumor <- function(spectrum,
                         sigs,
                         m.opts = NULL,
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
  
  # x0 is uniform distribution of mutations across signatures
  # (Each signature gets the same number of mutations)
  x0 <- rep(sum(spectrum) / ncol(sigs), ncol(sigs))
  
  if (any(is.na(x0))) {
    traceback()
    stop("Programming error, got an NA in x0 (global nloptr)")
  }
  global.res <- nloptr::nloptr(
    x0          = x0,
    eval_f      = m.opts[["global_eval_f"]],
    lb          = rep(0, ncol(sigs)),
    ub          = rep(sum(spectrum), ncol(sigs)),
    opts        = m.opts$global.opts,
    spectrum    = spectrum,
    sigs        = sigs,
    ...)
  if (m.opts$trace > 0) {
    message("global.res$objective = ", global.res$objective)
  }
  
  my.x0 <- global.res$solution
  
  # message("XX sum(spectrum) = ", sum(spectrum), "; sum(my.x0) = ", sum(my.x0))
  my.x0 <- my.x0 * (sum(spectrum) / sum(global.res$solution))
  # message("YY sum(spectrum) = ", sum(spectrum), "; sum(my.x0) = ", sum(my.x0))
  
  if (any(is.na(x0))) {
    traceback()
    stop("Programming error, got an NA in x0 (local nloptr)")
  }
  local.res <- nloptr::nloptr(
    x0          = my.x0,
    eval_f      = m.opts[["local_eval_f"]],
    eval_g_ineq = m.opts[["local_eval_g_ineq"]],
    lb          = rep(0, ncol(sigs)),
    ub          = rep(sum(spectrum) + 1e-2, ncol(sigs)),
    opts        = m.opts$local.opts,
    spectrum    = spectrum,
    sigs        = sigs,
    ...)
  if (m.opts$trace > 0)
    message("local.res$objective = ", local.res$objective)
  message("local.res$iterations = ", local.res$iterations)

  maxeval.warning <- NULL
  warnings <- NULL
  if (local.res$iterations == m.opts$local.opts[["maxeval"]]) {
    msg <- paste("reached maxeval on local optimization: ", local.res$iterations)
    if (m.opts$trace > 0) message(msg)
    warnings <- msg
  }

  names(local.res$solution) <- colnames(sigs)
  return(list(objective      =  local.res$objective,
              solution       = local.res$solution,
              warnings        = warnings,
              global.search.diagnostics = global.res,
              local.search.diagnostics = local.res))
}

