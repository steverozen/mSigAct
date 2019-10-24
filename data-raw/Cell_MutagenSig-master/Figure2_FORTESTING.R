dir.create("Figure2_results")
stopifnot(exists("SPIKE.TEST"))
setwd("./Figure2_results/")
source("../Header.R")
samples_details <- read.table("../00_data/final_mutagen_info_forR_v4_u.txt",sep = "\t",header = T,as.is = T, quote="\"")

sub_tab_all_info <- read.table("../00_data/denovo_subclone_subs_final.txt",sep = "\t",header = T, as.is = T)
sub_tab_all_info <- sub_tab_all_info[sub_tab_all_info$PM.Tum>=0.2,] # 172480/183133 = 0.94
sub_tab_all_info$Sample.Name <- sub("\\_.*","",sub_tab_all_info$Sample)
sub_tab_all_info <- sub_tab_all_info[sub_tab_all_info$Sample.Name!="MSM0",]
# SR sub_tab_all_info is a VCF-like file
stopifnot(dim(sub_tab_all_info) == c(172480, 18)) # SR

##########################
# Figure 2A 
# Substitutions
##########################

sub.summary <- data.frame(table(sub_tab_all_info$Sample))
# SR 2 column data frame mappig ids like MSM0.10_s1 to mutation counts, e.g. 278

names(sub.summary) <- c("Sample","sub_num")
sub.summary$Sample.Name <- sub("\\_.*","",sub.summary$Sample)
muts_summary_ddply <- ddply(sub.summary,c("Sample.Name"),summarise,NChild=length(sub_num),mean=mean(sub_num),sd=sd(sub_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name")

muts_summary_ddply_details_order <- muts_summary_ddply_details
#muts_summary_ddply_details_order[muts_summary_ddply_details_order$Group=="Control","Group"] <- "zControl"
muts_summary_ddply_details_order <- muts_summary_ddply_details_order[order(muts_summary_ddply_details_order$Group, muts_summary_ddply_details_order$mean, decreasing = T),]
#muts_summary_ddply_details_order$chemical_info <- paste0(muts_summary_ddply_details_order$Compound.Abbreviation," (", muts_summary_ddply_details_order$Concentration, ")")
#indel_mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#808080")


# Calculate pvalue
sub.summary <- data.frame(table(sub_tab_all_info$Sample)) # SR calculating it a second time?
names(sub.summary) <- c("Sample","sub_num")
sub.summary$Sample.Name <- sub("\\_.*","",sub.summary$Sample)
muts_summary_ddply <- ddply(sub.summary,c("Sample.Name"),summarise,NChild=length(sub_num),mean=mean(sub_num),sd=sd(sub_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name")
muts_summary_details <- merge(sub.summary, samples_details, by="Sample.Name")
write.table(muts_summary_details[,c(1,2,3,7,8)],"muts_summary_details.txt",sep = "\t",col.names = T, row.names = F, quote = F)

muts_summary_ddply_details_order <- muts_summary_ddply_details

muts_compound_ddply <- muts_summary_ddply_details_order[muts_summary_ddply_details_order$Group!="a_Control",]
Sample_withSig <- NULL

muts_control <- muts_summary_details[muts_summary_details$Group=="a_Control",]

################
# For compound have 2 subclones
# 
# SR Get distributon of means of samples of 2 controls.
sel_parents_mean_all_2 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],2,replace=T)  # SR dim(muts_control)[1] == nrow(muts_control)
  sel_parents <- muts_control[sel_idx,"sub_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_2 <- c(sel_parents_mean_all_2,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==2,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){ # SR for  i in 1:nrow(a), e.g for each compound.
  print(i)
  
  # SR get the empirical p value
  pvalue <- 
    round(
      length(which(sel_parents_mean_all_2 > round(a[i,"mean"]))) / 
        length(sel_parents_mean_all_2),
      digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all
#write.table(a,"subnum_distribution_2subclones.txt",sep = "\t",col.names = T, row.names = F, quote = F)

Sample_withSig <- a

################
# For compound have 3 subclones
# 
# SR Get distributon of means of samples of 3 controls.
sel_parents_mean_all_3 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],3,replace=T)
  sel_parents <- muts_control[sel_idx,"sub_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_3 <- c(sel_parents_mean_all_3,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==3,]

# ----- Replace MEG 1.25 uM with test data
meg.1.25.idx <- which(a$Sample.Name == "MSM0.124")
a[meg.1.25.idx, "mean"] <- mean(colSums(SPIKE.TEST$test.spectra))
a[meg.1.25.idx, "sd"] <- sd(colSums(SPIKE.TEST$test.spectra))

pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- 
    round(length(which(sel_parents_mean_all_3 > round(a[i,"mean"]))) /
            length(sel_parents_mean_all_3),
          digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all
#write.table(a,"subnum_distribution_3subclones.txt",sep = "\t",col.names = T, row.names = F, quote = F)
Sample_withSig <- rbind(Sample_withSig,a)

# For compound have 4 subclones
sel_parents_mean_all_4 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],4,replace=T)
  sel_parents <- muts_control[sel_idx,"sub_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_4 <- c(sel_parents_mean_all_4,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==4,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- round(length(which(sel_parents_mean_all_4>round(a[i,"mean"])))/length(sel_parents_mean_all_4),digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all
#write.table(a,"subnum_distribution_4subclones.txt",sep = "\t",col.names = T, row.names = F, quote = F)

Sample_withSig <- rbind(Sample_withSig,a)
Sample_withSig$P_adjust <- round(p.adjust(Sample_withSig$pvalue,method = "BH"),2)
Sample_withSig <- Sample_withSig[Sample_withSig$Sample.Name == "MSM0.124", ]


write.table(Sample_withSig,"Fig2A.tsv",sep = "\t",col.names = T, row.names = F, quote = F)

setwd("..")
# SR we could save key variables in an envirnment at this point

