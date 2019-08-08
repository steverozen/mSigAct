dir.create("Figure3_results")
setwd("./Figure3_results/")
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

###############################################
#
#    Sub catalogue (96 channels)
#
###############################################
sub_catalogue <- gen_muttype_new(sub_tab_all_info)
write.table(sub_catalogue,"sub_catalogue.txt",sep = "\t",col.names = T, row.names = F, quote = F)

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
for(j in 1:1000){
  chosen_control <- apply(control_profile,1,function(x) x[sample(1:length(x),1,replace=T)])
  
  chosen_control_all <- rbind(chosen_control_all,chosen_control)
}
chosen_control_muts <- t(chosen_control_all)
chosen_control_muts <- as.data.frame(chosen_control_muts)

centroid_control <- rowMeans(chosen_control_muts)


compoundlist <- data.frame(table(muts_compound$Sample.Name))
names(compoundlist) <- c("Sample.Name","Freq")
compoundlist$profile_SNR <- 0
for(i in 1:dim(compoundlist)[1]){
  
  print(i)
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
Sample_withSig_SNR <- read.table("../Figure2_results/Fig2A.txt",sep = "\t",header = T, as.is = T, quote="\"")
Sample_withSig_SNR2 <- merge(Sample_withSig_SNR,compoundlist[,c("Sample.Name","profile_SNR")],by="Sample.Name")

#############################################
# Extract signature for each compound condition
#############################################
control_profile <- sub_catalogue[,as.character(muts_control[,"Sample"])]
row.names(control_profile) <- sub_catalogue$MutationType

compoundlist <- data.frame(table(muts_compound$Sample.Name))
names(compoundlist) <- c("Sample.Name","Freq")
muttype_freq_template <- read.table("../00_data/MutationType_template.txt", sep = "\t", header = T, as.is = T)
stability_all <- NULL
#compoundlist$P_occurance <- 0
for(i in 1:dim(compoundlist)[1]){
  
  print(i)
  currentcompound <- sub_catalogue[,c(which(colnames(sub_catalogue) %in% muts_summary_details[muts_summary_details$Sample.Name==compoundlist[i,"Sample.Name"],"Sample"]))]
  row.names(currentcompound) <- sub_catalogue$MutationType
  
  currentcompound_stability <- ExtractSig_centroid_subs_noerrorbar(control_profile,currentcompound,1000,1.65,muttype_freq_template,4.2,7,paste0("subs_",as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Treatment"]))
  stability_all <- rbind(stability_all,c(currentcompound_stability))
  #MCSignature_subs_2_noerrorbar(control_profile,currentcompound,1000,1.65,4.2,7,paste0("subs_",as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Compound.Abbreviation"]))
  #compoundlist[i,"P_occurance"] <- current_P_occurance
  #ChemExposure(controlclones,compoundclones,6,5,paste0("indels_",as.character(compoundlist[i,1]),"_",currentcompound[1,"Compound.Abbreviation"],"_exposure"))
}
stability_all_2 <- cbind(stability_all,compoundlist)
stability_all_2 <- as.data.frame(stability_all_2)
names(stability_all_2) <- c("min_simi","max_simi","Sample.Name","subclone_num")
Sample_withSig_SNR3 <- merge(Sample_withSig_SNR2,stability_all_2,by="Sample.Name")
write.table(Sample_withSig_SNR3,"Fig3B.txt",sep = "\t",col.names = T, row.names = F, quote = F)

#############################################
# heatmap Fig3C
#############################################
Sample_withSig <- read.table("./Fig3B.txt",sep = "\t",header = T, as.is = T,quote = "\"")
sig_subs_full <- Sample_withSig[Sample_withSig$pvalue<=0.01 & Sample_withSig$profile_SNR>=2,]
sig_all <- NULL
for(i in 1:dim(sig_subs_full)[1]){
  sig_file <- read.table(paste0("subs_",sig_subs_full[i,"Sample.Name"],"_",sig_subs_full[i,"Treatment"],"_exposure.txt"), sep = "\t", header = T, as.is = T)
  sig_all <- cbind(sig_all,sig_file[,"percentage"])
  
}
control_sig <- read.table("./Fig3A.txt", sep = "\t", header = T, as.is = T)
control_sig$percentage <- control_sig$mean/sum(control_sig$mean)
sig_all <- cbind(sig_all,control_sig$percentage)
sig_all <- t(as.data.frame(sig_all))
row.names(sig_all) <- c(paste0(sig_subs_full$Compound.Abbreviation,"_",sig_subs_full$Concentration),"control")

colnames(sig_all) <- sig_file$MutationType
sig_all <- data.matrix(sig_all)
sig_all <- round(sig_all,2)
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "red"))(n = 128)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(0,  # for red
               0.05,           # for yellow
               0.2,max(sig_all))             # for green

pdf(file="Fig3C.pdf", h=10, w=15, onefile=TRUE)
gplots::heatmap.2(sig_all,
                  hclustfun=function(x) hclust(x,method = 'complete'),
                  reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), # Reorder dendrogram by branch means rather than sums
                  main = "Heatmap of substitution signature", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat map
                  margins =c(12,9),     # widens margins around plot
                  keysize=1,
                  col=my_palette,       # use on color palette defined earlier
                  key = TRUE, 
                  #      breaks=col_breaks,    # enable color transition at specified limits
                  dendrogram="row",     # only draw a row dendrogram
                  Colv="NA")            # turn off column clustering
dev.off()



