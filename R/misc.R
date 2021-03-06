# Needed to silence devtools::check warnings about gobal variable used in dplyr
utils::globalVariables(c(".data"))


List2TibbleRow <- function(a.list) {

  list2 <- lapply(a.list, function(x) {
    if (length(x) > 1) {
      return(list(x))
    } else if (length(x) == 0) {
      return(NA)
    } else {
      return(x)
    }})

  rr <- do.call(tibble::tibble, list2)
  return(rr)
}

ListOfList2Tibble <- function(list.of.lists) {
  rr <- lapply(list.of.lists, List2TibbleRow)
  do.call("rbind", rr)
}


if (FALSE) {
  testin <- list(list(a = 1, b = 2), list(a = 3, b = 33))
}

#' Given signatures (sigs) and exposures (exp), return a spectrum or spectra
#'
#' @param sigs Signature as a matrix or or data frame, with each
#'      row one mutation type (g.e. CCT > CAT or CC > TT) and
#'      each column a signature.
#'
#' @param exp The exposures for one or more samples as a matrix or data.frame,
#'      with each row a signature and each column a sample.
#'
#' @param use.sig.names If \code{TRUE} check that \code{rownames(exp)} is
#'   a subset of \code{colnames(sigs)}, and use only the columns in \code{sigs}
#'   that are present in \code{exp}.
#'
# @return The matrix product \code{sigs %*% exp} after some error checking.
#'
#' @details Does not care or check if \code{colSums(sigs) == 1}.
#'   Error checking is minimal since this function is called often.
#'
#' @export
#'
ReconstructSpectrum <- function(sigs, exp, use.sig.names = FALSE) {
  stopifnot(is.matrix(sigs))
  if (!(is.matrix(exp))) {
     exp <- as.matrix(exp)
  }
  if (use.sig.names) {
    sigs <- sigs[ , rownames(exp)]
  }
  stopifnot(length(exp) == ncol(sigs))
  return(as.matrix(sigs) %*% as.matrix(exp))
}


#' For backward compatibility.
#'
#' Use \code{\link{ReconstructSpectrum}} instead.
#'
#' @inheritParams ReconstructSpectrum
#'
# @return The matrix product \code{sigs %*% exp}.
#'
#' @keywords internal

prop.reconstruct <- function(sigs, exp) {
  return(ReconstructSpectrum(sigs, exp, FALSE))
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


#' Cosine similarity with useful argument types..
#'
#' Calls \code{\link[lsa]{cosine}}.
#'
#' @param v1 A vector or single-column matrix
#' @param v2 A vector or single-column matrix
#'
#' @export

cossim <- function(v1, v2) {
  if (!is.null(ncol(v1)))  {
    stopifnot(ncol(v1) == 1)
    v1 <- v1[ , 1]
  }
  if (!is.null(ncol(v2)))  {
    stopifnot(ncol(v2) == 1)
    v2 <- v2[ , 1]
  }
  lsa::cosine(v1, v2)

}

ClosestCosSig <- function(spectrum) {
  cos <-
    apply(PCAWG7::signature$genome$SBS96,
          MARGIN = 2,
          FUN =
            function(sig) {
              lsa::cosine(as.vector(sig), as.vector(spectrum))})
  max.cos <- which(cos == max(cos))
  return(cos[max.cos])

}


ClosestCosSigDensity <- function(spectrum) {
  spec <-
    ICAMS::TransformCatalog(
      spectrum,
      target.catalog.type = "density")

  sigs <-
    ICAMS::TransformCatalog(
      PCAWG7::signature$genome$SBS96,
      target.catalog.type = "density.signature")

  cos <-
    apply(sigs,
          MARGIN = 2,
          FUN =
            function(sig) {
              lsa::cosine(as.vector(sig), as.vector(spec))})
  max.cos <- which(cos == max(cos))
  return(cos[max.cos])

}


LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}


ClearWarnings <- function() {
  assign("last.warning", NULL, envir = baseenv())
}


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


#' Euclidean reconstruction error.
#'
#' @keywords internal
EDist2SpectRounded <- function(exp, sig.names, spect) {
  reconstruction <-
    prop.reconstruct(
      sigs = PCAWG7::signature$genome$SBS96[ , sig.names],
      exp = round(exp))
  reconstruction <- round(reconstruction)
  class(spect) <- "matrix"
  err <- stats::dist(t(cbind(reconstruction, spect)), method = "euclidean")
  return(err)
}


#' Euclidean reconstruction error.
#'
#' @keywords internal
EDist2Spect <- function(exp, sig.names, spect) {
  reconstruction <-
    prop.reconstruct(
      sigs =
        PCAWG7::signature$genome$SBS96[ , sig.names], exp = exp)
  class(spect) <- "matrix"
  err <- stats::dist(t(cbind(reconstruction, spect)), method = "euclidean")
  return(err)
}


Adj.mc.cores <- function(mc.cores) {
  if (Sys.info()["sysname"] == "Windows" && mc.cores > 1) {
    message("On Windows, changing mc.cores from ", mc.cores, " to 1")
    return(1)
  }
  return(mc.cores)
}


