#!Rscript

# Script to create a 1-sample spectrum catalog from the PCAWG7 
# data package, from the PCAWG platinum data.

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
  stop("Usage: get.spect.R <mutation type> <sample id>")
}

mutation.type <- args[1]
id <- args[2]
message(mutation.type)

known.sig.types <-c("DBS78","ID", "SBS96", "SBS192", "SBS1536")
if (!mutation.type %in% known.sig.types) {
  stop("First argument must be one of ",
          paste(known.sig.types, collapse = ", "))
}

message(id)

# "Bladder-TCC::SP1003"

ss <- PCAWG7::spectra$PCAWG[[mutation.type]]
if (!id %in% colnames(ss)) stop(id, " not found")

cc <- ss[ , id, drop = FALSE]

id <- gsub("::", "_", id)
ff <- paste0(mutation.type, "_", id, ".csv")

cc <- ICAMS::as.catalog(cc)

ICAMS::WriteCatalog(cc, file = ff)