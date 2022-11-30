#' @keywords internal
NullReturnForMAPAssignActivity <- 
  function(null.assignment, null.spect, msg, time.used = system.time(3)) {
    return(
      list(proposed.assignment           = null.assignment,
           proposed.reconstruction       = null.spect,
           reconstruction.distances      = NULL,
           all.tested                    = NULL,
           alt.solutions                 = NULL,
           time.for.MAP.assign           = time.used,
           error.messages                = msg
      ))
  }