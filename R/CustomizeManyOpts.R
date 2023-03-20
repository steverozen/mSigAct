CustomizeManyOpts <- function(loglh.fn) {
  object.fn <- function(exp, spectrum, sigs) {
    reconstruction <- sigs %*% exp
    loglh <- loglh.fn(spectrum = spectrum, expected.counts = reconstruction)
    return(-loglh)
  }

  g_ineq_for_object.fn <- function(exp, spectrum, sigs) {
    retval <- abs(sum(exp) - sum(spectrum))
    return(retval)
  }

  return(list(
    global.opts       = DefaultGlobalOpts(),
    local.opts        = DefaultLocalOpts(),
    trace             = 0,
    global_eval_f     = object.fn,
    local_eval_f      = object.fn,
    local_eval_g_ineq = g_ineq_for_object.fn
  ))
}
