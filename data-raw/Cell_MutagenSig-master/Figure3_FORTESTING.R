if (!dir.exists("Figure3_results")) dir.create("Figure3_results")
setwd("./Figure3_results/")


###############################################
#
#    Sub catalogue (96 channels)
#
###############################################
###############################################


#----- SPIKE IN TEST DATA

stopifnot(dim(sub_tab_all_info) == c(172480, 18)) # SR
#  sub_catalogue <- gen_muttype_new(sub_tab_all_info) # CAUTION row order is
#  "A[C>A]A" "A[C>A]C" "A[C>A]G" "A[C>A]T" "A[C>G]A" ,
#  NOT
#  A[C>A]A A[C>A]C A[C>A]G A[C>A]T C[C>A]A (mSigAct::kucab.sigs, see
#  kucab.sigs <- read.table("../../kucab.sigs.from.mendely.2019.10.25/Mutagen53_sub_signature.txt", sep = "\t", header = T)
#  

# kucab.sub.catalog <- sub_catalogue
# usethis::use_data(K_sub_catalogue)
sub_catalogue <- mSigAct::kucab.sub.catalog
stopifnot(dim(sub_catalogue) == c(96, 326)) # SR

# write.table(sub_catalogue,"sub_catalogue.txt",sep = "\t",col.names = T, row.names = F, quote = F)

meg.idx2 <- grep("MSM0.124", colnames(sub_catalogue)) # SR
sub_catalogue[  , meg.idx2] <- SPIKE.TEST$test.spectra
write.table(sub_catalogue,"sub_catalogue_spike.tsv",sep = "\t",col.names = T, row.names = F, quote = F)

#----- END SPIKE IN TEST DATA

muts_control <- muts_summary_details[muts_summary_details$Group=="a_Control",]
muts_compound <- muts_summary_details[muts_summary_details$Group!="a_Control",]

# Averaged control profile
parentmuts <- sub_catalogue[,c(1,which(colnames(sub_catalogue) %in% muts_summary_details[muts_summary_details$Group=="a_Control","Sample"]))]
control_mean <- plotCountbasis_average_se(parentmuts,2,7,"Fig3A")

# Mutational signatures
#############################################
#   Calculate signal-to-noise ratio distance for mutational profile
#    between mutagen and distinction (control)
#############################################
control_profile <- sub_catalogue[,as.character(muts_control[,"Sample"])]

chosen_control_all <- NULL
for(j in 1:1000) {
  chosen_control <- apply(control_profile,1,function(x) x[sample(1:length(x),1,replace=T)])
  
  chosen_control_all <- rbind(chosen_control_all,chosen_control)
}
chosen_control_muts <- t(chosen_control_all)
chosen_control_muts <- as.data.frame(chosen_control_muts)

centroid_control <- rowMeans(chosen_control_muts)


compoundlist <- data.frame(table(muts_compound$Sample.Name))
names(compoundlist) <- c("Sample.Name","Freq")
compoundlist <- compoundlist[compoundlist$Sample.Name == "MSM0.124", ]
compoundlist$profile_SNR <- 0
for(i in 1:dim(compoundlist)[1]){
  
  c.name <- compoundlist[i, "Sample.Name"]
  cat(i, c.name, "\n")
  if (c.name != "MSM0.124") { 
    cat("skipping")
    next
  }
  
  currentcompound <- sub_catalogue[,as.character(muts_compound[muts_compound$Sample.Name==as.character(compoundlist[i,"Sample.Name"]),"Sample"])]
  
  # permutation
  chosen_compound_all <- NULL
  for(j in 1:1000){
    chosen_compound <- apply(currentcompound,1,function(x) x[sample(1:length(x),1,replace=T)])
    chosen_compound_all <- rbind(chosen_compound_all,chosen_compound)
  }
  
  chosen_compound_all <- t(chosen_compound_all)
  chosen_compound_all <- as.data.frame(chosen_compound_all)
  
  
  centroid_compound <- rowMeans(chosen_compound_all)
  sd_compound <- sd_highD(chosen_compound_all,centroid_compound-centroid_control)
  sd_control<- sd_highD(chosen_control_muts,centroid_compound-centroid_control)
  
  compoundlist[i,"profile_SNR"]  <- norm(as.matrix(centroid_compound-centroid_control),"f")/(sd_compound+sd_control)
}

Sample_withSig_SNR <- Sample_withSig

Sample_withSig_SNR2 <- merge(Sample_withSig_SNR,compoundlist[,c("Sample.Name","profile_SNR")],by="Sample.Name")

#############################################
# Extract signature for each compound condition
#############################################
control_profile <- sub_catalogue[,as.character(muts_control[,"Sample"])]
row.names(control_profile) <- sub_catalogue$MutationType

# compoundlist <- data.frame(table(muts_compound$Sample.Name))
# names(compoundlist) <- c("Sample.Name","Freq")
muttype_freq_template <- read.table("../00_data/MutationType_template.txt", sep = "\t", header = T, as.is = T)
stability_all <- NULL
#compoundlist$P_occurance <- 0
for(i in 1:dim(compoundlist)[1]){
  
  
  c.name <- compoundlist[i, "Sample.Name"]
  cat(i, c.name,"\n")
  if (c.name != "MSM0.124") { 
    cat("skipping")
    next
  }
  
  
  
  currentcompound <- sub_catalogue[,c(which(colnames(sub_catalogue) %in% muts_summary_details[muts_summary_details$Sample.Name==compoundlist[i,"Sample.Name"],"Sample"]))]
  row.names(currentcompound) <- sub_catalogue$MutationType
  
  currentcompound_stability <- ExtractSig_centroid_subs_noerrorbar(control_profile,currentcompound,1000,1.65,muttype_freq_template,4.2,7,paste0("subs_",as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Treatment"]))
  stability_all <- rbind(stability_all,c(currentcompound_stability))
}
stability_all_2 <- cbind(stability_all,compoundlist)
stability_all_2 <- as.data.frame(stability_all_2)
names(stability_all_2) <- c("min_simi","max_simi","Sample.Name","subclone_num")
Sample_withSig_SNR3 <- merge(Sample_withSig_SNR2,stability_all_2,by="Sample.Name")
write.table(Sample_withSig_SNR3,"Fig3B-test.txt",sep = "\t",col.names = T, row.names = F, quote = F)
print(Sample_withSig_SNR3)
test.sig.table <- 
  read.table("subs_MSM0.124_Methyleugenol (1.25 mM)_exposure.txt",
             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(lsa::cosine(test.sig.table$percentage, me.sig.ICAMS), "\n")
# 0.9633224  ratio = 0.5 add.noise = FALSE
# 0.9570034  ratio = 0.5 add.noise = TRUE
# 0.6619694   ratio = 0.25 add.noise = TRUE, fails on p_adjust, fails on SNR, fails on max_simi
# test.sig.table rows are in the *STANDARD* order
setwd("..")
