# Spectra from HepG2 exposed to various nitrosamines

options(error = browser)
options(warn = 3)
# showMethods("getSeq")
# Function: getSeq (package Biostrings)
# x="BSgenome"
# x="FaFile"
# x="FaFileList"
# x="TwoBitFile"
# x="XStringSet"

# Important
SBS.dir <- devtools::package_file(
  file.path("data-raw", "nitro-vcfs-2020-01-09"))
tmp.files <- grep("\\.vcf$",
                 list.files(SBS.dir,
                            full.names = TRUE), 
                 ignore.case = TRUE, value = TRUE) 
vcf.names <- sub(".*(N..._cl.).*", "\\1", tmp.files)

nitrosamine.examples <- 
  ICAMS::StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
    files = tmp.files,
    ref.genome = 
      BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    trans.ranges = ICAMS::trans.ranges.GRCh37,
    region = "genome",
    names.of.VCFs = vcf.names,
    output.file = file.path(SBS.dir, "nitrosamines")
  )

usethis::use_data(nitrosamine.examples, overwrite = TRUE)
rm(tmp.files, SBS.dir, vcf.names, nitrosamine.examples)


