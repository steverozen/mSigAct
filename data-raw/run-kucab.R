

setwd("data-raw/Cell_MutagenSig-master/")
# source("Figure3_SubstitutionSig.R")
# 
new.dir <- paste0("Fig3_", gsub("[ :]", "_", date()))
dir.create(new.dir)
# Copied directly from Figure3_SubstitionsSig.R:
setwd(new.dir)
source("../Header.R")
samples_details <- read.table("../00_data/final_mutagen_info_forR_v4_u.txt",sep = "\t",header = T,as.is = T, quote="\"")

sub_tab_all_info <- read.table("../00_data/denovo_subclone_subs_final.txt",sep = "\t",header = T, as.is = T)
sub_tab_all_info <- sub_tab_all_info[sub_tab_all_info$PM.Tum>=0.2,] # 172480/183133 = 0.94
sub_tab_all_info$Sample.Name <- sub("\\_.*","",sub_tab_all_info$Sample)
sub_tab_all_info <- sub_tab_all_info[sub_tab_all_info$Sample.Name!="MSM0",]

sub.summary <- data.frame(table(sub_tab_all_info$Sample))
names(sub.summary) <- c("Sample","sub_num")
sub.summary$Sample.Name <- sub("\\_.*","",sub.summary$Sample)
muts_summary_ddply <- ddply(sub.summary,c("Sample.Name"),summarise,NChild=length(sub_num),mean=mean(sub_num),sd=sd(sub_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name")
muts_summary_details <- merge(sub.summary, samples_details, by="Sample.Name")

sub_catalogue <- gen_muttype_new(sub_tab_all_info)
# Commented out in mSigActwrite.table(sub_catalogue,"sub_catalogue.txt",sep = "\t",col.names = T, row.names = F, quote = F)

# New code
kucab.spectra <- sub_catalogue
rownames(kucab.spectra) <- ICAMS:::Unstaple96(kucab.spectra$MutationType)
kucab.spectra <- kucab.spectra[ICAMS::catalog.row.order$SBS96, ]
kucab.spectra <- kucab.spectra[ , -(1:2)]
rownames(muts_summary_details) <- muts_summary_details$Sample
new.col.names.table <- muts_summary_details[colnames(kucab.spectra),
                                            c("Treatment", "S9.mix", "Group", "Sample")]

new.col.names <- apply(new.col.names.table, MARGIN = 1, FUN = paste, collapse = "_")
new.col.names <- sub(" ", "_", new.col.names)
colnames(kucab.spectra) <- new.col.names
muts_summary_details$New.names <- new.col.names

kucab.spectra <- 
  ICAMS::as.catalog(kucab.spectra, 
                    region = "genome", 
                    catalog.type = "counts",
                    ref.genome = "hg19")

ICAMS::PlotCatalogToPdf(kucab.spectra[ , sort(colnames(kucab.spectra))], "kucab.spectra.pdf")

# control.index <- which(muts_summary_details$Group == "a_Control")

muts_control <- muts_summary_details[muts_summary_details$Group == "a_Control",]
muts_compound <- muts_summary_details[muts_summary_details$Group != "a_Control",]

# START HERE
# Averaged control profile
parentmuts <- kucab.spectra[ ,muts_control$New.names]
pmut <- cbind(MutationType = rownames(parentmuts), parentmuts)

# Throws error:
# control_mean <- plotCountbasis_average_se(pmut,2,7,"Fig3A")

# Mutational signatures
#############################################
#   Calculate signal-to-noise ratio distance for mutational profile
#    between mutagen and distinction (control)
#############################################
# control_profile <- sub_catalogue[,as.character(muts_control[,"Sample"])]
# End of material copied from Kucab


# Step 1, get the control spectra from the Kucab et al. data.
# The variable control_profile contains the control spectra
# and is defined in Figure3_SubstitutionSig.R.
kucab.controls <- parentmuts

# rownames(kucab.controls) <- ICAMS:::Unstaple96(rownames(kucab.controls))
# kucab.controls <- 
#   ICAMS::as.catalog(kucab.controls, catalog.type = "counts", region = "genome")
# control.details <- 
#   muts_summary_details[muts_summary_details$Group == "a_Control", ]
# new.colnames <- 
#    cbind(colnames(kucab.controls),
#     control.details[colnames(kucab.controls),
#                     c("Compound", "Concentration")])
# colnames(kucab.controls) <- 
#  apply(new.colnames, MARGIN = 1, FUN = paste0, collapse = "_")
# kucab.controls <- kucab.controls[ICAMS::catalog.row.order$SBS96, ]
# 
# usethis::use_data(kucab.controls, overwrite = TRUE)
ICAMS::PlotCatalogToPdf(kucab.controls, "kucab.controls.pdf")


# Step 2, compute the distribution of means of samples of
# spectra from 2, 3, and 4 cell lines.
kucab.compute.control.dist <- function(num.controls) {
  x <- colSums(mSigAct::kucab.controls)
  retval <- NULL
  set.seed(217)
  for(j in 1:10000){
    selx <- sample(x, num.controls, replace = TRUE)
    retval <- c(retval, mean(selx))
  }
  return(retval)
}

kucab.control.dist <- list()
for (i in 2:4) kucab.control.dist[[i]] <- kucab.compute.control.dist(i)
usethis::use_data(kucab.control.dist)





if (FALSE) {
# The functions below are temporary -- have been moved to ICAMS
XUnstaple96 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           substr(c1, 3, 3),
           substr(c1, 7, 7),
           substr(c1, 5, 5))
  return(retval)
}

XRestaple96 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           "[",
           substr(c1, 2, 2),
           ">",
           substr(c1, 4, 4),
           "]",
           substr(c1, 3, 3))
  return(retval)
}
}
