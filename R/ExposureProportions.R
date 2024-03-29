#' Return the proportions of tumors of a given cancer type that have a particular signature
#'
#' @param mutation.type A character string, one of "SBS96", "SBS192", "ID", 
#'    "DBS78".
#'
#' @param cancer.type A character string. For some common cancer types, see
#'   \link[PCAWG7]{CancerTypes} for more details.
#'
#' @param all.sigs An optional matrix of known signatures, 
#'    with column names being signatures ids. Only used to drop 
#'    signatures not present in \code{all.sigs}.
#'
#' @param drop.sigs.no.info If TRUE, drop signatures not present in 
#'    the column names of \code{all.sigs}.
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
#' @examples 
#' cancer.types <- PCAWG7::CancerTypes()
#' cancer.types
#' sigs.prop <- ExposureProportions(mutation.type = "SBS96", 
#'                                  cancer.type = "Lung-AdenoCA")
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
      sig.names <- PCAWG7::SBS96_ID_to_SBS192_ID(sig.names)
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
