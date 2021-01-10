#' Return the proportions of tumors of a given cancer type that have a particular signature
#'
#' @param mutation.type A character string, one of "SBS96", "SBS192", "ID", "DBS78".
#'
#' @param cancer.type A character string.
#'
#' @param all.sigs An optional matrix of known signatures, 
#'    with column names being signatures ids. Only used to drop 
#'    signatures not present in \code{all.sigs}.
#'
#' @param drop.sigs.no.info If TRUE, drop signatures any not present in 
#'    the column names of \code{all.sigs}.
#'    There are some signatures that do not have SBS192 versions, including SBS29.
#'    
#' @param must.include A character vector of signature IDs that
#'    must be included, even if they have not previously been
#'    observed in that cancer type. The associated proportion is
#'    specified by \code{must.include.prop}.
#'    
#' @param must.include.prop The value used for the expected proportion of
#'    signatures in \code{must.include} but not previously observed
#'    in the given \code{cancer.type}.
#'
#' @export
#'
#' @return A numerical vector of the proportion of 
#'   tumors of type \code{cancer.type} with each signature
#'   for those signatures observed in \code{cancer.type}. 
#'   The names are the signature ids.
#'
ExposureProportions <- function(
  mutation.type, 
  cancer.type, 
  all.sigs          = NULL, 
  drop.sigs.no.info = TRUE,
  must.include      = character(),
  must.include.prop  = 0.1
  ) {

  if (cancer.type %in% c("Unknown", "None")) {
    sigs.prop <- numeric()
    sig.names <- character()
  } else {
    # We don't have exposure information for SBS192 signatures,
    # so we use the exposure information for the SBS96 signatures
    # with the same names.
    mtype <- ifelse(mutation.type == "SBS192", "SBS96", mutation.type)
    sigs.prop <- PCAWG7::exposure.stats$PCAWG[[mtype]][[cancer.type]]
    if (is.null(sigs.prop)) {
      stop("Cannot find exposure information for ",
           mtype, " for ", cancer.type)
    }
    rm(mtype)
    sig.names <- rownames(sigs.prop)
    sigs.prop <- unlist(sigs.prop[ , 2])
    if (mutation.type == "SBS192") {
      sig.names <- AddMinusEFor192(sig.names)
    }
    names(sigs.prop) <- sig.names
  }
  
  other.prop.sig.names <- setdiff(must.include, sig.names)
  sig.names <- c(sig.names, other.prop.sig.names)
  other.props <-
    rep.int(must.include.prop, times = length(other.prop.sig.names))
  names(other.props) <- other.prop.sig.names
  sigs.prop <- c(sigs.prop, other.props)
  
  if (!is.null(all.sigs)) {
    sigs.no.info <- setdiff(sig.names, colnames(all.sigs))

    for (zz in sigs.no.info) {
      message("ExposureProportions:")
      message("No signature ", zz, " for ", mutation.type, 
              " in all.sigs\n(",
              paste(colnames(all.sigs), collapse = ", "),
              ")")
      message("Dropping from sigs.prop")
      sigs.prop <- sigs.prop[-(which(names(sigs.prop) == zz))]
    }

    if (!drop.sigs.no.info) 
      stop("There were signaures with no proporiton information")
  }
  
  return(sigs.prop)
}

#' Translate SBS96 signature IDs to SBS192 signature IDs.
#'
#' @param sig.ids Character vector of SBS96 signature IDs.
#' 
#' @keywords internal
#' @return Character vector of corresponding SBS192 signature IDs; some
#'   have "-E" (for exome) post-pended.
AddMinusEFor192 <- function(sig.ids) {
  need.minus.E <- 
    c("SBS23", "SBS25", "SBS27", "SBS29",
      "SBS42", "SBS45", "SBS46", "SBS59")
  rr <- lapply(sig.ids, 
               function(x) ifelse(x %in% need.minus.E, paste0(x, "-E"), x))
  return(unlist(rr))
}

#' Return a character vector of the IDs of possible SBS96 signature artifacts.
#'
#' @export
PossibleArtifacts <- function() {
  return(c("SBS27", "SBS43", "SBS45", "SBS46", "SBS47",
           "SBS48", "SBS49", "SBS50", "SBS51", "SBS52",
           "SBS53", "SBS54", "SBS55", "SBS56",
           "SBS57", "SBS58", "SBS59", "SBS60"))
}

#' Return a character vector of the IDs of rare SBS96 signatures.
#'
#' @export
RareSignatures <- function() {
  return(paste0("SBS", 84:90))
}
