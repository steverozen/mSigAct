#' @keywords internal
check.mclapply.result <- function(ss, where, names = NULL) {
  ok <- TRUE
  for (i in 1:length(ss)) {
    
    if (is.null(names)) {
      name <- paste("item", i)
    } else {
      name = names[i]
    }
    if (is.null(ss[[i]])) {
      message(where, ": got NULL return for ", name)
      ok <- FALSE
    }
    if ("try-error" %in% class(ss[[i]])) {
      message(where, ": got try-error return for ", name)
      print(ss[i])
      ok <- FALSE
    }
  }
  stopifnot(ok)
}
