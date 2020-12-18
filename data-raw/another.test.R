library(BSgenome.Hsapiens.1000genomes.hs37d5)

spect <- 
  ICAMS::MutectVCFFilesToCatalog(
    "data-raw/mSigAct-example-VCFs/mSigAct-example-VCFs/Mutect-mSigAct-example-GRCh37-s1.vcf",
    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5)


# incomplete


rr <- OneMAPAssignTest(
  spect = spect$catSBS96,
  reference.exp = PCAWG7::spectra
)