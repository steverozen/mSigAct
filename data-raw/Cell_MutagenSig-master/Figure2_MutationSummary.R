dir.create("Figure2_results")
setwd("./Figure2_results/")
source("../Header.R")
samples_details <- read.table("../00_data/final_mutagen_info_forR_v4_u.txt",sep = "\t",header = T,as.is = T, quote="\"")

sub_tab_all_info <- read.table("../00_data/denovo_subclone_subs_final.txt",sep = "\t",header = T, as.is = T)
sub_tab_all_info <- sub_tab_all_info[sub_tab_all_info$PM.Tum>=0.2,] # 172480/183133 = 0.94
sub_tab_all_info$Sample.Name <- sub("\\_.*","",sub_tab_all_info$Sample)
sub_tab_all_info <- sub_tab_all_info[sub_tab_all_info$Sample.Name!="MSM0",]

##########################
# Figure 2A 
# Substitutions
##########################

sub.summary <- data.frame(table(sub_tab_all_info$Sample))
names(sub.summary) <- c("Sample","sub_num")
sub.summary$Sample.Name <- sub("\\_.*","",sub.summary$Sample)
muts_summary_ddply <- ddply(sub.summary,c("Sample.Name"),summarise,NChild=length(sub_num),mean=mean(sub_num),sd=sd(sub_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name")

muts_summary_ddply_details_order <- muts_summary_ddply_details
#muts_summary_ddply_details_order[muts_summary_ddply_details_order$Group=="Control","Group"] <- "zControl"
muts_summary_ddply_details_order <- muts_summary_ddply_details_order[order(muts_summary_ddply_details_order$Group, muts_summary_ddply_details_order$mean, decreasing = T),]
#muts_summary_ddply_details_order$chemical_info <- paste0(muts_summary_ddply_details_order$Compound.Abbreviation," (", muts_summary_ddply_details_order$Concentration, ")")
#indel_mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#808080")
mypalette <- c("#808080","#fabebe","#46f0f0","#aa6e28","#f032e6","#008080","#d2f53c","#e6beff","#911eb4","#ffe119","#0082c8","#f58231","#3cb44b","#e6194b")


pdf(file="Fig2A.pdf", onefile=TRUE,width=15,height=5)
p <- ggplot(data=muts_summary_ddply_details_order, aes(x=Sample.Name, y=mean, fill=Group))+ geom_bar(stat="identity",position="dodge")+xlab("Sample")+ylab("Count")
p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",position=position_dodge(.9),width=.1)+scale_y_continuous()
p <- p+scale_x_discrete(limits = muts_summary_ddply_details_order$Sample.Name,labels = muts_summary_ddply_details_order$Treatment)
p <- p+scale_fill_manual(values=mypalette)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black",hjust = 1),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),
             axis.title.y = element_text(size=15),
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))

print(p)
dev.off()
#write.table(muts_summary_ddply_details_order, paste0("sub_summary_ddply_details",".txt"),sep = "\t",col.names = T, row.names = F, quote = F)


# Calculate pvalue
sub.summary <- data.frame(table(sub_tab_all_info$Sample))
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
# For compound have 2 subclones
sel_parents_mean_all_2 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],2,replace=T)
  sel_parents <- muts_control[sel_idx,"sub_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_2 <- c(sel_parents_mean_all_2,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==2,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- round(length(which(sel_parents_mean_all_2>round(a[i,"mean"])))/length(sel_parents_mean_all_2),digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all
#write.table(a,"subnum_distribution_2subclones.txt",sep = "\t",col.names = T, row.names = F, quote = F)

Sample_withSig <- a
# For compound have 3 subclones
sel_parents_mean_all_3 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],3,replace=T)
  sel_parents <- muts_control[sel_idx,"sub_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_3 <- c(sel_parents_mean_all_3,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==3,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- round(length(which(sel_parents_mean_all_3>round(a[i,"mean"])))/length(sel_parents_mean_all_3),digits = 4)
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

write.table(Sample_withSig,"Fig2A.txt",sep = "\t",col.names = T, row.names = F, quote = F)




##########################
# Figure 2B
# Double Substitutions
##########################
DinucleotidesList <- FindDinucleotides(sub_tab_all_info)
samplelist <- data.frame(table(sub_tab_all_info$Sample))
names(samplelist) <- c("Sample","sub_num")
dinuc <- data.frame(table(DinucleotidesList$Sample))
names(dinuc) <- c("Sample","dinuc_num")
sample_dinuc_list <- merge(samplelist,dinuc,by="Sample",all.x=T)
sample_dinuc_list[is.na(sample_dinuc_list)] <- 0
sample_dinuc_list$dinuc_percentage <- sample_dinuc_list$dinuc_num/sample_dinuc_list$sub_num
sample_dinuc_list$Sample.Name <- sub("\\_.*","",sample_dinuc_list$Sample)
sample_dinuc_list <- merge(sample_dinuc_list, samples_details, by="Sample.Name")
write.table(sample_dinuc_list,"Dinucleotides_number_summary.txt",sep = "\t",col.names = T, row.names = F, quote = F)
sample_dinuc_list <- read.table("Dinucleotides_number_summary.txt", sep = "\t", header = T, as.is = T,quote = "\"")
sample_dinuc_list$adjust_dinuc_num <- sample_dinuc_list$dinuc_num/2
muts_summary_ddply <- ddply(sample_dinuc_list,c("Sample.Name"),summarise,NChild=length(adjust_dinuc_num),mean=mean(adjust_dinuc_num),sd=sd(adjust_dinuc_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name")

muts_summary_ddply_details_order <- muts_summary_ddply_details
#muts_summary_ddply_details_order[muts_summary_ddply_details_order$Group=="Control","Group"] <- "zControl"
muts_summary_ddply_details_order <- muts_summary_ddply_details_order[order(muts_summary_ddply_details_order$Group, muts_summary_ddply_details_order$mean, decreasing = T),]
#muts_summary_ddply_details_order$chemical_info <- paste0(muts_summary_ddply_details_order$Compound.Abbreviation," (", muts_summary_ddply_details_order$Concentration, ")")
# grey     pink   light blue   brown      marg  dark green   lime   light purple  purple  yellow     blue      orange   green      red 
mypalette <- c("#808080","#fabebe","#46f0f0","#aa6e28","#f032e6","#008080","#d2f53c","#e6beff","#911eb4","#ffe119","#0082c8","#f58231","#3cb44b","#e6194b")
filename <- paste0("Fig2B", ".pdf")
pdf(file=filename, onefile=TRUE,width=15,height=5)
p <- ggplot(data=muts_summary_ddply_details_order, aes(x=Sample.Name, y=mean, fill=Group))+ geom_bar(stat="identity",position="dodge")+xlab("Sample")+ylab("Count")
p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",position=position_dodge(.9),width=.1)+scale_y_continuous()
p <- p+scale_x_discrete(limits = muts_summary_ddply_details_order$Sample.Name,labels = muts_summary_ddply_details_order$Treatment)
p <- p+scale_fill_manual(values=mypalette)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black",hjust=1),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),
             axis.title.y = element_text(size=15),
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))

print(p)
dev.off()

# Calculate pvalue
Dinucleotides_profile <- DinucleotidesList[DinucleotidesList$neigbor_dist==-1,]
Dinucleotides_profile$dinuc_Ref <- paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Ref_neighbor)
Dinucleotides_profile$dinuc_Alt <- paste0(Dinucleotides_profile$Alt,Dinucleotides_profile$Alt_neighbor)
Dinucleotides_profile$dinuc_mutation <- paste0(Dinucleotides_profile$dinuc_Ref,">",Dinucleotides_profile$dinuc_Alt)

Dinucleotides_profile$dinuc_Ref_rc <- as.character(reverseComplement(DNAStringSet(Dinucleotides_profile$dinuc_Ref)))
Dinucleotides_profile$dinuc_Alt_rc <- as.character(reverseComplement(DNAStringSet(Dinucleotides_profile$dinuc_Alt)))

Dinucleotides_profile$dinuc_Ref_final <- ifelse(Dinucleotides_profile$dinuc_Ref<=Dinucleotides_profile$dinuc_Ref_rc,Dinucleotides_profile$dinuc_Ref,Dinucleotides_profile$dinuc_Ref_rc)
Dinucleotides_profile$dinuc_Alt_final <- ifelse(Dinucleotides_profile$dinuc_Ref<=Dinucleotides_profile$dinuc_Ref_rc,Dinucleotides_profile$dinuc_Alt,Dinucleotides_profile$dinuc_Alt_rc)
Dinucleotides_profile$dinuc_mutation_final <- paste0(Dinucleotides_profile$dinuc_Ref_final,">",Dinucleotides_profile$dinuc_Alt_final)

sub.summary <- data.frame(table(Dinucleotides_profile$Sample))
names(sub.summary) <- c("Sample","dinu_num")

subs_summary_allsample <- data.frame(table(sub_tab_all_info$Sample))
names(subs_summary_allsample) <- c("Sample","subs")
muts_summary_details <- merge(sub.summary, subs_summary_allsample, by="Sample", all=T)
muts_summary_details[is.na(muts_summary_details)] <- 0
muts_summary_details$Sample.Name <- sub("\\_.*","",muts_summary_details$Sample)
muts_summary_details <- merge(muts_summary_details, samples_details, by="Sample.Name", all=T)

muts_summary_details <- muts_summary_details[order(muts_summary_details$Group, muts_summary_details$Sample.Name),]
write.table(muts_summary_details,"dinu_summary_details.txt",sep = "\t",col.names = T, row.names = F, quote = F)

muts_summary_ddply <- ddply(muts_summary_details,c("Sample.Name"),summarise,NChild=length(dinu_num),mean=mean(dinu_num),sd=sd(dinu_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name", all=T)
muts_summary_ddply_details[is.na(muts_summary_ddply_details)] <- 0


muts_summary_ddply_details_order <- muts_summary_ddply_details

muts_compound_ddply <- muts_summary_ddply_details_order[muts_summary_ddply_details_order$Group!="a_Control",]
Sample_withSig <- NULL

muts_control <- muts_summary_details[muts_summary_details$Group=="a_Control",]

# For compound have 2 subclones
sel_parents_mean_all_2 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],2,replace=T)
  sel_parents <- muts_control[sel_idx,"dinu_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_2 <- c(sel_parents_mean_all_2,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==2,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- round(length(which(sel_parents_mean_all_2>round(a[i,"mean"])))/length(sel_parents_mean_all_2),digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all

Sample_withSig <- a
# For compound have 3 subclones
sel_parents_mean_all_3 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],3,replace=T)
  sel_parents <- muts_control[sel_idx,"dinu_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_3 <- c(sel_parents_mean_all_3,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==3,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- round(length(which(sel_parents_mean_all_3>round(a[i,"mean"])))/length(sel_parents_mean_all_3),digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all
Sample_withSig <- rbind(Sample_withSig,a)

# For compound have 4 subclones
sel_parents_mean_all_4 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],4,replace=T)
  sel_parents <- muts_control[sel_idx,"dinu_num"]
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

Sample_withSig <- rbind(Sample_withSig,a)
Sample_withSig$P_adjust <- p.adjust(Sample_withSig$pvalue,method = "BH")
write.table(Sample_withSig,"Fig2B.txt",sep = "\t",col.names = T, row.names = F, quote = F) # for double nucleotides, significant: the mean > 3 and qvalue<=0.01

##########################
# Figure 2C
# Indels
##########################
indel.classified <- read.table("../00_data/denovo_subclone_indels.final.txt",sep = "\t",header = T, as.is = T)
indel.classified$Sample.Name <- sub("\\_.*","",indel.classified$Sample)
indel.classified <- indel.classified[indel.classified$Sample.Name!="MSM0",]
indel.classified_details <- merge(indel.classified, samples_details, by="Sample.Name")
indel.classified_details <- indel.classified_details[indel.classified_details$VAF.Tum_Cal>=0.2,]

indel.summary <- data.frame(table(indel.classified_details$Sample))
names(indel.summary) <- c("Sample","indel_num")
indel.summary$Sample.Name <- sub("\\_.*","",indel.summary$Sample)
muts_summary_ddply <- ddply(indel.summary,c("Sample.Name"),summarise,NChild=length(indel_num),mean=mean(indel_num),sd=sd(indel_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name")

muts_summary_ddply_details_order <- muts_summary_ddply_details
muts_summary_ddply_details_order <- muts_summary_ddply_details_order[order(muts_summary_ddply_details_order$Group, muts_summary_ddply_details_order$mean, decreasing = T),]
write.table(muts_summary_ddply_details_order, paste0("indel_summary_ddply_details",".txt"),sep = "\t",col.names = T, row.names = F, quote = F)

               # grey     pink   light blue   brown      marg  dark green   lime   light purple  purple  yellow     blue      orange   green      red 
mypalette <- c("#808080","#fabebe","#46f0f0","#aa6e28","#f032e6","#008080","#d2f53c","#e6beff","#911eb4","#ffe119","#0082c8","#f58231","#3cb44b","#e6194b")
filename <- paste0("Fig2C", ".pdf")
pdf(file=filename, onefile=TRUE,width=15,height=5)
p <- ggplot(data=muts_summary_ddply_details_order, aes(x=Sample.Name, y=mean, fill=Group))+ geom_bar(stat="identity",position="dodge")+xlab("Sample")+ylab("Count")
p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",position=position_dodge(.9),width=.1)+scale_y_continuous(limits=c(0,210),breaks=(seq(0,200,40)))
p <- p+scale_x_discrete(limits = muts_summary_ddply_details_order$Sample.Name,labels = (muts_summary_ddply_details_order$Treatment))
p <- p+scale_fill_manual(values=mypalette)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black",hjust = 1),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),
             axis.title.y = element_text(size=15),
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))

print(p)
dev.off()


# Calculate p value
indel.summary <- data.frame(table(indel.classified_details$Sample))
names(indel.summary) <- c("Sample","indel_num")
indel.summary$Sample.Name <- sub("\\_.*","",indel.summary$Sample)
muts_summary_ddply <- ddply(indel.summary,c("Sample.Name"),summarise,NChild=length(indel_num),mean=mean(indel_num),sd=sd(indel_num),se=sd/sqrt(NChild))
muts_summary_ddply_details <- merge(muts_summary_ddply, samples_details, by="Sample.Name")

muts_summary_ddply_details_order <- muts_summary_ddply_details
muts_summary_details <- merge(indel.summary, samples_details, by="Sample.Name")

muts_compound_ddply <- muts_summary_ddply_details_order[muts_summary_ddply_details_order$Group!="a_Control",]
Sample_withSig <- NULL

# control
muts_control <- muts_summary_details[muts_summary_details$Group=="a_Control",]

# For compound have 2 subclones
sel_parents_mean_all_2 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],2,replace=T)
  sel_parents <- muts_control[sel_idx,"indel_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_2 <- c(sel_parents_mean_all_2,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==2,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- round(length(which(sel_parents_mean_all_2>round(a[i,"mean"])))/length(sel_parents_mean_all_2),digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all

Sample_withSig <- a
# For compound have 3 subclones
sel_parents_mean_all_3 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],3,replace=T)
  sel_parents <- muts_control[sel_idx,"indel_num"]
  sel_parents_mean <- mean(sel_parents)
  sel_parents_mean_all_3 <- c(sel_parents_mean_all_3,sel_parents_mean)
}
a <- muts_compound_ddply[muts_compound_ddply$NChild==3,]
pvalue_all <- NULL
for(i in 1:dim(a)[1]){
  print(i)
  pvalue <- round(length(which(sel_parents_mean_all_3>round(a[i,"mean"])))/length(sel_parents_mean_all_3),digits = 4)
  pvalue_all <- c(pvalue_all,pvalue)
}
a$pvalue <- pvalue_all
Sample_withSig <- rbind(Sample_withSig,a)

# For compound have 4 subclones
sel_parents_mean_all_4 <- NULL
for(j in 1:10000){
  sel_idx <- sample(1:dim(muts_control)[1],4,replace=T)
  sel_parents <- muts_control[sel_idx,"indel_num"]
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

Sample_withSig <- rbind(Sample_withSig,a)
Sample_withSig$P_adjust <- round(p.adjust(Sample_withSig$pvalue,method = "BH"),2)

write.table(Sample_withSig,"Fig2C.txt",sep = "\t",col.names = T, row.names = F, quote = F)


##############################################################################
# Figure 2D
# Mutagenic potential (Mutagenicity ratio between control and treatment)
##############################################################################

# subs
sub.summary <- data.frame(table(sub_tab_all_info$Sample))
names(sub.summary) <- c("Sample","sub_num")
sub.summary$Sample.Name <- sub("\\_.*","",sub.summary$Sample)
sub.summary_details <- merge(sub.summary, samples_details, by="Sample.Name")

control_num <- sub.summary_details[sub.summary_details$Group=="a_Control","sub_num"]

Sample_withSig <- read.table("./Fig2A.txt",sep = "\t",header = T, as.is = T,quote = "\"")
sig_subs_full <- Sample_withSig[Sample_withSig$P_adjust<=0.01,]
sig_all <- NULL
for(i in 1:dim(sig_subs_full)[1]){
  mutagen_num <- sub.summary_details[sub.summary_details$Sample.Name==sig_subs_full[i,"Sample.Name"],"sub_num"]
  sig_all <- rbind(sig_all,Mutagenicity_ratio(control_num,mutagen_num))
}
sig_all <- as.data.frame(sig_all)
names(sig_all) <- c("ratio","ratio_se")
sig_all$Sample.Name <- sig_subs_full$Sample.Name
sig_all_details <- merge(sig_all, samples_details, by="Sample.Name")
sig_all_details <- sig_all_details[order(sig_all_details$ratio, decreasing=T),]
sig_all_details$treatment <- paste0(sig_all_details$Sample.Name,"_",sig_all_details$Compound.Abbreviation,"_",sig_all_details$Concentration)
write.table(sig_all_details,"Fig2D_subs.txt",sep = "\t",col.names = T, row.names = F, quote = F)



# Double substitutions
muts_summary_details <- read.table("dinu_summary_details.txt",sep = "\t",header = T, as.is = T,quote = "\"")
control_num <- muts_summary_details[muts_summary_details$Group=="a_Control","dinu_num"]

Sample_withSig <- read.table("Fig2B.txt",sep = "\t",header = T, as.is = T,quote = "\"")
sig_indel_full <- Sample_withSig[Sample_withSig$P_adjust<=0.01 & Sample_withSig$mean>3,]

sig_all <- NULL
for(i in 1:dim(sig_indel_full)[1]){
  mutagen_num <- muts_summary_details[muts_summary_details$Sample.Name==sig_indel_full[i,"Sample.Name"],"dinu_num"]
  #sig_all <- rbind(sig_all,(mean(mutagen_num)-mean(control_num))/mean(control_num))
  sig_all <- rbind(sig_all,Mutagenicity_ratio(mean(control_num),mutagen_num))
  
}
sig_all <- as.data.frame(sig_all)
names(sig_all) <- c("ratio","ratio_se")
sig_all$Sample.Name <- sig_indel_full$Sample.Name
sig_all_details <- merge(sig_all, samples_details, by="Sample.Name")
sig_all_details <- sig_all_details[order(sig_all_details$ratio, decreasing=T),]
write.table(sig_all_details,"Fig2D_doublesubs.txt",sep = "\t",col.names = T, row.names = F, quote = F)



# indels
indel.summary <- data.frame(table(indel.classified_details$Sample))
names(indel.summary) <- c("Sample","indel_num")
indel.summary$Sample.Name <- sub("\\_.*","",indel.summary$Sample)
indel.summary_details <- merge(indel.summary, samples_details, by="Sample.Name")

control_num <- indel.summary_details[indel.summary_details$Group=="a_Control","indel_num"]
Sample_withSig <- read.table("Fig2C.txt",sep = "\t",header = T, as.is = T,quote = "\"")
sig_indel_full <- Sample_withSig[Sample_withSig$P_adjust<=0.01,]

sig_all <- NULL
for(i in 1:dim(sig_indel_full)[1]){
  mutagen_num <- indel.summary_details[indel.summary_details$Sample.Name==sig_indel_full[i,"Sample.Name"],"indel_num"]
  sig_all <- rbind(sig_all,Mutagenicity_ratio(control_num,mutagen_num))
}
sig_all <- as.data.frame(sig_all)
names(sig_all) <- c("ratio","ratio_se")
sig_all$Sample.Name <- sig_indel_full$Sample.Name
sig_all_details <- merge(sig_all, samples_details, by="Sample.Name")
sig_all_details <- sig_all_details[order(sig_all_details$ratio, decreasing=T),]
write.table(sig_all_details,"Fig2D_indels.txt",sep = "\t",col.names = T, row.names = F, quote = F)

# Combine subs, dinus, indels together
sig_all_details_subs <- read.table("Fig2D_subs.txt",sep = "\t",header = T, as.is = T)
sig_all_details_subs$type <- "1sub"
sig_all_details_dinu <- read.table("Fig2D_doublesubs.txt",sep = "\t",header = T, as.is = T)
sig_all_details_dinu$type <- "2dinu"
sig_all_details_2 <- rbind(sig_all_details_subs[,c("Sample.Name","ratio","ratio_se","type")],sig_all_details_dinu[,c("Sample.Name","ratio","ratio_se","type")])
sig_all_details_indels <- read.table("Fig2D_indels.txt",sep = "\t",header = T, as.is = T)
sig_all_details_indels$type <- "3indel"
sig_all_details_2 <- rbind(sig_all_details_2,sig_all_details_indels[,c("Sample.Name","ratio","ratio_se","type")])
a <- data.frame(table(sig_all_details_2$Sample.Name))
names(a) <- c("Sample.Name","freq")
sig_all_details_2 <- merge(sig_all_details_2,a,by="Sample.Name")
sig_all_details_2 <- merge(sig_all_details_2,samples_details,by="Sample.Name")


sig_all_details_2[sig_all_details_2$freq==1 & sig_all_details_2$type=="1sub","freq"] <- "sub"
sig_all_details_2[sig_all_details_2$freq==1 & sig_all_details_2$type=="3indel","freq"] <- "indel"
sig_all_details_2[sig_all_details_2$freq==2,"freq"] <- "sub_indel"
sig_all_details_2[sig_all_details_2$freq==3,"freq"] <- "all"

sig_all_details_2 <- sig_all_details_2[order(sig_all_details_2$Group,sig_all_details_2$Compound.Abbreviation,sig_all_details_2$type),]
treatment_order <- unique(sig_all_details_2$Treatment)
filename=paste0("Fig2D",".pdf")
pdf(file=filename, onefile=TRUE,width = 6,height = 20)
g1 <-ggplot(sig_all_details_2, aes(x=Treatment, y=type)) + geom_tile(aes(fill=ratio),colour="white")+geom_text(aes(label=paste(round(ratio,digits = 1))),size=5)
g1 <-g1 +scale_fill_gradient2(high="red", low="white",space="Lab")
#g1 <-g1 +scale_x_discrete(limits = unique(as.character(cossimil$sample1)))
#g1 <-g1 +scale_y_discrete(limits = unique(as.character(cossimil$sample2[order(as.character(cossimil$sample2))])))
g1 <- g1+scale_x_discrete(limits = treatment_order)
g1 <-g1 +theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10,colour = "black",hjust=1),
               axis.text.y=element_text(size=10,colour = "black"),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
g1 <- g1+coord_flip()
print(g1)
dev.off()

