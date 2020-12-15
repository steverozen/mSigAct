#!Rscript
args <- commandArgs(trailingOnly = TRUE)

mutation.type <- args[1]
id <- args[2]
message(mutation.type)

if (!mutation.type %in% c("DBS78"))

message(id)

# "Bladder-TCC::SP1003"

ss <- PCAWG7::spectra$PCAWG[[mutation.type]]
if (!id %in% colnames(ss)) stop(id, " not found")

cc <- ss[ , id, drop = FALSE]

id <- gsub("::", "_", id)
ff <- paste0(mutation.type, "_", id, ".csv")

cc <- ICAMS::as.catalog(cc)

ICAMS::WriteCatalog(cc, file = ff)