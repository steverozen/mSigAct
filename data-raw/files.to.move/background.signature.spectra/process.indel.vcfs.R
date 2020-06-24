# Create background indel spectra catalogs and plot them

working.dir <- 
  devtools::package_file(
    "data-raw/background.signature.spectra/backgroundVCFS-2020-02-04/")

vcf.names <- 
  list.files(
    path = working.dir,
    pattern = "*._INDEL.*",
    full.names = TRUE)

short.names <- sub("_INDEL.*", "", sub(".*/", "", vcf.names))

cats <- ICAMS::StrelkaIDVCFFilesToCatalog(
  files = vcf.names,
  ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
  region = "genome",
  names.of.VCFs = short.names)

ICAMS::WriteCatalog(cats$catalog, file = file.path(working.dir, "background.indels.csv"))

ICAMS::PlotCatalogToPdf(cats$catalog, file = file.path(working.dir, "background.indels.pdf"))

                        