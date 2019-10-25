# plot kucab indels

library(data.table)
pre.vcf <- fread("00_data/denovo_subclone_indels.final.txt")
nm <- toupper(colnames(pre.vcf))
colnames(pre.vcf) <- nm
vcf.list <- split(x = pre.vcf, f = pre.vcf$SAMPLE)
names(vcf.list)
vcf18 <- vcf.list[c("MSM0.18_s1", "MSM0.18_s2")]
all.id.cat <- ICAMS::VCFsToIDCatalogs(vcf.list, ref.genome = "hg19", region = "genome")
ICAMS::PlotCatalogToPdf(all.id.cat, "all.id.cat.pdf")
