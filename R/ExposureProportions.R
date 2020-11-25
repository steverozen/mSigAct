#' Return the proportions of tumors of a given cancer type that have a particular signature
#'
#' @param mutation.type A character string, one of "SBS96", "SBS192", "ID", "DBS78".
#'
#' @param cancer.type A character string.
#'
#' @param all.sigs An optional matrix of known signatures, with column names being signatures ids.
#'
#' @param drop.sigs.no.info If TRUE, drop any not present in the column names of \code{all.sigs}.
#'    There are some signatures that do not have SBS192 versions, including SBS29.
#'
#' @export
#'
#' @return A numerical vector of the proportion of tumors of type \code{cancer.type} with each signature
#'   for those signatures observed in \code{cancer.types}.The names are the signature ids.
#'
ExposureProportions <- function(mutation.type, cancer.type, all.sigs = NULL, drop.sigs.no.info = TRUE) {

  mtype <- ifelse(mutation.type == "SBS192", "SBS96", mutation.type)

  sigs.prop <- PCAWG7::exposure.stats$PCAWG[[mtype]][[cancer.type]]
  if (is.null(sigs.prop)) {
    stop("Cannot find exposure information for ",
         mtype, " for ", cancer.type)
  }
  sig.names <- rownames(sigs.prop)
  sigs.prop <- unlist(sigs.prop[ , 2])
  names(sigs.prop) <- sig.names

  if (!is.null(all.sigs)) {
    sigs.no.info <- setdiff(sig.names, colnames(all.sigs))

    for (zz in sigs.no.info) {
      message("No signature ", zz, " for ", mutation.type)
      message("Dropping from signature universe")
      sigs.prop <- sigs.prop[-(which(names(sigs.prop) == zz))]
    }

    if (!drop.sigs.no.info) stop("There were signaures with no proporiton information")
  }
  return(sigs.prop)
}

