dir.create("Figure4_results")
setwd("./Figure4_results/")
source("../Header.R")
samples_details <- read.table("../00_data/final_mutagen_info_forR_v4_u.txt",sep = "\t",header = T,as.is = T, quote="\"")


#############
# Figure 4A
#############
sub_tab_all_info <- read.table("../00_data/denovo_subclone_subs_final.txt",sep = "\t",header = T, as.is = T)
sub_tab_all_info <- sub_tab_all_info[sub_tab_all_info$PM.Tum>=0.2,] # 172480/183133 = 0.94

DinucleotidesList <- FindDinucleotides(sub_tab_all_info)
write.table(DinucleotidesList,"DinucleotidesList.txt",sep = "\t",col.names = T, row.names = F, quote = F)

samplelist <- data.frame(table(sub_tab_all_info$Sample))
names(samplelist) <- c("Sample","sub_num")
dinuc <- data.frame(table(DinucleotidesList$Sample))
names(dinuc) <- c("Sample","dinuc_num")
sample_dinuc_list <- merge(samplelist,dinuc,by="Sample",all.x=T)
sample_dinuc_list[is.na(sample_dinuc_list)] <- 0
sample_dinuc_list$dinuc_percentage <- sample_dinuc_list$dinuc_num/sample_dinuc_list$sub_num
sample_dinuc_list$Sample.Name <- sub("\\_.*","",sample_dinuc_list$Sample)
sample_dinuc_list <- merge(sample_dinuc_list, samples_details, by="Sample.Name")
DinucleotidesList <- sample_dinuc_list
# The probability of no dinucleiotides, taking context into consideration


  TrinucleotideFreq_origin <- read.table("../00_data/trinucleotides_counts_genome.txt", sep = "\t", header = T, as.is = T)
  TrinucleotideFreq_origin$avail_freq <- TrinucleotideFreq_origin$freq
  TrinucleotideFreq <- TrinucleotideFreq_origin
  
  subs_context2bp <- read.table("../00_data/subs.txt",sep = "\t", header = T, as.is = T)
  subs_context2bp <- subs_context2bp[,c("VariantID","pre_context","rear_context")]
  names(subs_context2bp) <- c("VariantID","pre_context2","rear_context2")
  sub_tab_all_info2 <- merge(sub_tab_all_info, subs_context2bp, by="VariantID")
  samplelist <- data.frame(table(sub_tab_all_info2$Sample))
  names(samplelist) <- c("Sample","sub_num")
  samplelist$P_nodinuc <- 0
  for(i in 1:dim(samplelist)[1]){
    print(i)
    subs <- sub_tab_all_info2[sub_tab_all_info2$Sample==samplelist[i,"Sample"],]
    P_nodinuc <- 1
    TrinucleotideFreq <- TrinucleotideFreq_origin
    
    for(j in 1:dim(subs)[1]){
      
      sub_context <- paste0(subs[j,"pre_context"], subs[j,"Ref"], subs[j,"rear_context"])
      P_nodinuc <- P_nodinuc*TrinucleotideFreq[TrinucleotideFreq$trinuc==sub_context,"avail_freq"]/TrinucleotideFreq[TrinucleotideFreq$trinuc==sub_context,"freq"]
      # adjust TrinucleotideFreq
      TrinucleotideFreq[TrinucleotideFreq$trinuc==sub_context,"freq"] <- TrinucleotideFreq[TrinucleotideFreq$trinuc==sub_context,"freq"] -1
      TrinucleotideFreq[TrinucleotideFreq$trinuc==sub_context,"avail_freq"] <- TrinucleotideFreq[TrinucleotideFreq$trinuc==sub_context,"avail_freq"] -1
      
      neigbor_context_1 <- paste0(subs[j,"Ref"], subs[j,"rear_context2"])
      neigbor_context_2 <- paste0(subs[j,"pre_context2"], subs[j,"Ref"])
      TrinucleotideFreq[TrinucleotideFreq$trinuc==neigbor_context_1,"avail_freq"] <- TrinucleotideFreq[TrinucleotideFreq$trinuc==neigbor_context_1,"avail_freq"] -1
      TrinucleotideFreq[TrinucleotideFreq$trinuc==neigbor_context_2,"avail_freq"] <- TrinucleotideFreq[TrinucleotideFreq$trinuc==neigbor_context_2,"avail_freq"] -1
    }
    samplelist[i,"P_nodinuc"] <- P_nodinuc
  }
  samplelist$P_dinuc <- 1-samplelist$P_nodinuc
  samplelist$P_dinuc_adjust <- p.adjust(samplelist$P_dinuc,method = "BH")
  
  DinucleotidesList_P <- merge(DinucleotidesList,samplelist[,c("Sample","P_nodinuc","P_dinuc","P_dinuc_adjust")], by="Sample")
  write.table(DinucleotidesList_P,"Fig4A.txt",sep = "\t",col.names = T, row.names = F, quote = F)
  
  DinucleotidesList_P <- read.table("Fig4A.txt", sep = "\t", header = T, as.is = T,quote = "\"")
  DinucleotidesList_P$adjust_dinuc_num <- DinucleotidesList_P$dinuc_num/2
  DinucleotidesList_P$Group <- factor(DinucleotidesList_P$Group, levels = rev(levels(factor(DinucleotidesList_P$Group))))
                        
  indel_mypalette <- c("#e6194b", # red
                       "#3cb44b", # green
                       "#f58231", # orange
                       "#0082c8", # blue
                       "#ffe119", # yellow
                       "#911eb4", # purple
                       "#e6beff", # light purple
                       "#d2f53c", # light green
                       "#008080", # dark green
                       "#f032e6", # magenta
                       "#aa6e28", # brown
                       "#46f0f0", # light blue
                       "#fabebe", # pink
                       "#808080") # grey
  filename <- paste0("Fig4A", ".pdf")
  pdf(file=filename, onefile=TRUE,width=6.6,height=4,useDingbats=FALSE)
  p <- ggplot(data=DinucleotidesList_P, aes(x=sub_num, y=P_dinuc_adjust))+ geom_point(aes(colour=Group))+xlab("Number of subs")+ylab("Adjusted pvalue")
  # p <- p+scale_y_continuous(breaks=seq(0, 0.01, 0.005))
  p <- p+scale_colour_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
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




#############
# Figure 4B
#############
  # Dinucleotides profile
DinucleotidesList <- read.table("DinucleotidesList.txt", sep = "\t", header = T, as.is = T)
Dinucleotides_profile <- DinucleotidesList[DinucleotidesList$neigbor_dist==-1,]
Dinucleotides_profile$dinuc_Ref <- paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Ref_neighbor)
Dinucleotides_profile$dinuc_Alt <- paste0(Dinucleotides_profile$Alt,Dinucleotides_profile$Alt_neighbor)
Dinucleotides_profile$dinuc_mutation <- paste0(Dinucleotides_profile$dinuc_Ref,">",Dinucleotides_profile$dinuc_Alt)
write.table(Dinucleotides_profile,"Dinucleotides.txt",sep = "\t",col.names = T, row.names = F, quote = F)

# heatmap
dimut <- c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G")

dinuc_MutationType <- NULL
dinuc_firstmut <- NULL
dinuc_second <- NULL
for(i in 1:length(dimut)){
  for(j in 1:length(dimut)){
    dinuc_MutationType <- c(dinuc_MutationType,paste0(dimut[i],":",dimut[j]))
    dinuc_firstmut <- c(dinuc_firstmut,dimut[i])
    dinuc_second <- c(dinuc_second,dimut[j])
  }
}
dinuc_template <- cbind(dinuc_MutationType,dinuc_firstmut)
dinuc_template <- data.frame(cbind(dinuc_template,dinuc_second))
names(dinuc_template) <- c("MutationType","dinuc_firstmut","dinuc_secondmut")

Dinucleotides_profile <- read.table("Dinucleotides.txt", sep = "\t", header = T, as.is = T)
a <- data.frame(table(Dinucleotides_profile$Sample.Name))
Dinucleotides_profile <- Dinucleotides_profile[Dinucleotides_profile$Sample.Name %in% a[a$Freq>=20,1],]


Dinucleotides_profile$firstmut <- ifelse(paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Alt)<=paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))),paste0(Dinucleotides_profile$Ref,">",Dinucleotides_profile$Alt),paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),">",as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))))
Dinucleotides_profile$secondmut <- ifelse(paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Alt)<=paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))),paste0(Dinucleotides_profile$Ref_neighbor,">",Dinucleotides_profile$Alt_neighbor),paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref))),">",as.character(complement(DNAStringSet(Dinucleotides_profile$Alt)))))
Dinucleotides_profile$dimut <- paste0(Dinucleotides_profile$firstmut, ":",Dinucleotides_profile$secondmut)
sigfile_freq <- data.frame(table(Dinucleotides_profile$Sample.Name,Dinucleotides_profile$dimut))
names(sigfile_freq) <- c("Sample.Name","MutationType","Freq")
control_sigset <-dcast(sigfile_freq,MutationType~Sample.Name,value.var = "Freq")
control_sigset <- merge(dinuc_template,control_sigset,by="MutationType",all.x=T)
control_sigset[is.na(control_sigset)] <- 0
sigfile_freq <- melt(control_sigset,id=c("MutationType","dinuc_firstmut","dinuc_secondmut"))
sigfile_freq$group <- 0
sigfile_freq[sigfile_freq$value>0,]$group <- 2^(floor(log2(sigfile_freq[sigfile_freq$value>0,"value"])))
names(sigfile_freq) <- c("MutationType","dinuc_firstmut","dinuc_secondmut","Sample.Name","Freq","Group")
sigfile_freq <- merge(sigfile_freq,samples_details[,c("Sample.Name","Compound.Abbreviation","Concentration")],by="Sample.Name",all.x=T)
sigfile_freq$treatment <- paste0(sigfile_freq$Sample.Name, "_",sigfile_freq$Compound.Abbreviation,"_",sigfile_freq$Concentration)
pdf(file=paste0("Fig4B_treatment_dinuc_profile_with20Dinuc.pdf"), onefile=TRUE,width=12,height=6)
p <- ggplot(data=sigfile_freq, aes(x=dinuc_firstmut, y=dinuc_secondmut,fill=factor(Group)))+ geom_tile()+xlab("firstmut")+ylab("secondmut")+ggtitle("mutagen_dinuc_profile")+scale_fill_brewer(palette="Purples")#+scale_fill_manual(values=c("#FFFFFF", "#CCCC00","#FFCC33","#FF9966","#FF6633","#CC3300","#660000"))
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=9,colour = "black"),
             axis.text.y=element_text(size=9,colour = "black"),
             axis.title.x = element_text(size=15),
             axis.title.y = element_text(size=15),
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))
p <- p+facet_wrap(~treatment,ncol=4, scale="free")
print(p)
dev.off()

# For dinucleotide number > 20 mutagens, Calculate cosimilarity
Dinucleotides_profile <- read.table("Dinucleotides.txt", sep = "\t", header = T, as.is = T)
a <- data.frame(table(Dinucleotides_profile$Sample.Name))
Dinucleotides_profile <- Dinucleotides_profile[Dinucleotides_profile$Sample.Name %in% a[a$Freq>=20,1],]

Dinucleotides_profile$firstmut <- ifelse(paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Alt)<=paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))),paste0(Dinucleotides_profile$Ref,">",Dinucleotides_profile$Alt),paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),">",as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))))
Dinucleotides_profile$secondmut <- ifelse(paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Alt)<=paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))),paste0(Dinucleotides_profile$Ref_neighbor,">",Dinucleotides_profile$Alt_neighbor),paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref))),">",as.character(complement(DNAStringSet(Dinucleotides_profile$Alt)))))
Dinucleotides_profile$dimut <- paste0(Dinucleotides_profile$firstmut, ":",Dinucleotides_profile$secondmut)
sigfile_freq <- data.frame(table(Dinucleotides_profile$Sample.Name,Dinucleotides_profile$dimut))
names(sigfile_freq) <- c("Sample.Name","MutationType","Freq")
sigfile_freq <- merge(sigfile_freq,samples_details[,c("Sample.Name","Compound.Abbreviation","Concentration")],by="Sample.Name",all.x=T)
sigfile_freq$treatment <- paste0(sigfile_freq$Sample.Name, "_",sigfile_freq$Compound.Abbreviation,"_",sigfile_freq$Concentration)

control_sigset <-dcast(sigfile_freq,MutationType~treatment,value.var = "Freq")
control_sigset <- merge(dinuc_template,control_sigset,by="MutationType",all.x=T)
control_sigset[is.na(control_sigset)] <- 0
control_sigset_seletcted <- control_sigset

cossimil <- as.matrix(proxy::simil(t(control_sigset_seletcted[,4:dim(control_sigset_seletcted)[2]]), diag=FALSE, upper=FALSE, method="cosine",auto_convert_data_frames = FALSE))
cossimil[lower.tri(cossimil,diag=TRUE)]=NA
cossimil <- as.data.frame(as.table(cossimil))
cossimil=na.omit(cossimil)

write.table(cossimil,"Fig4B_cossim_treatmentswith20Dinu.txt",sep="\t",col.names = T, row.names = F, quote = F)
names(cossimil) <- c("sample1","sample2","simil")
cossimil$simil <- round(cossimil$simil,2)
filename=paste0("Fig4B_cossim_treatmentswith20Dinu",".pdf")
pdf(file=filename, onefile=TRUE,width = 7.5,height = 6)
g1 <-ggplot(cossimil, aes(x=sample1, y=sample2)) + geom_tile(aes(fill=simil),colour="white")+geom_text(aes(label=paste(simil)),size=3)
g1 <-g1 +scale_fill_gradient2(midpoint=0.5,high="#B2182B", low="#2166AC",space="Lab",limits=c(0, 1))
g1 <-g1 +theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10),
               axis.text.y=element_text(size=10,colour = "black"),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
print(g1)
dev.off()

# Hcluster

dinuc_with20Dinuc <- control_sigset_seletcted[,4:dim(control_sigset_seletcted)[2]]
dinuc_with20Dinuc_t <- t(dinuc_with20Dinuc)
row.names(dinuc_with20Dinuc_t) <- colnames(dinuc_with20Dinuc)

dinuc_with20Dinuc_t_per <- dinuc_with20Dinuc_t/rowSums(dinuc_with20Dinuc_t)[row(dinuc_with20Dinuc_t)]
clusters <- hclust(stats::dist(dinuc_with20Dinuc_t_per, method = "euclidean"), method = "complete")
pdf(file="Fig4B_Dinuc_hclust_wide_20dinu.pdf", h=10, w=25, onefile=TRUE)
plot(clusters, hang = -1, cex=1.1)
dev.off()

  
  
#############
# Figure 4C
#############
# treatments have more than 20 dinucleotides,78 channels
Dinucleotides_profile <- read.table("./Dinucleotides.txt", sep = "\t", header = T, as.is = T)
Dinucleotides_profile$dinuc_Ref_rc <- as.character(reverseComplement(DNAStringSet(Dinucleotides_profile$dinuc_Ref)))
Dinucleotides_profile$dinuc_Alt_rc <- as.character(reverseComplement(DNAStringSet(Dinucleotides_profile$dinuc_Alt)))
  
Dinucleotides_profile$dinuc_Ref_final <- ifelse(Dinucleotides_profile$dinuc_Ref<=Dinucleotides_profile$dinuc_Ref_rc,Dinucleotides_profile$dinuc_Ref,Dinucleotides_profile$dinuc_Ref_rc)
Dinucleotides_profile$dinuc_Alt_final <- ifelse(Dinucleotides_profile$dinuc_Ref<=Dinucleotides_profile$dinuc_Ref_rc,Dinucleotides_profile$dinuc_Alt,Dinucleotides_profile$dinuc_Alt_rc)
  
Dinucleotides_profile$dinuc_mutation_final <- paste0(Dinucleotides_profile$dinuc_Ref_final,">",Dinucleotides_profile$dinuc_Alt_final)
  
a <- data.frame(table(Dinucleotides_profile$Sample.Name))
Dinucleotides_profile_20 <- Dinucleotides_profile[Dinucleotides_profile$Sample.Name %in% a[a$Freq>=20,1],]
sample_dinu <- gen_muttype_Dinucleotides(Dinucleotides_profile_20,"Sample.Name")
write.table(sample_dinu,"Fig4C_doublesubs_8signature.txt",sep = "\t",col.names = T, row.names = F, quote = F)
plot_allsample_Dinucleotides_profile(sample_dinu,2,13,20,"Fig4C.pdf")

#############
# Figure 4D
#############
Dinucleotides_profile <- read.table("Dinucleotides.txt", sep = "\t", header = T, as.is = T)
Dinucleotides_profile$dinuc_Ref_rc <- as.character(reverseComplement(DNAStringSet(Dinucleotides_profile$dinuc_Ref)))
Dinucleotides_profile$dinuc_Alt_rc <- as.character(reverseComplement(DNAStringSet(Dinucleotides_profile$dinuc_Alt)))

Dinucleotides_profile$dinuc_Ref_final <- ifelse(Dinucleotides_profile$dinuc_Ref<=Dinucleotides_profile$dinuc_Ref_rc,Dinucleotides_profile$dinuc_Ref,Dinucleotides_profile$dinuc_Ref_rc)
Dinucleotides_profile$dinuc_Alt_final <- ifelse(Dinucleotides_profile$dinuc_Ref<=Dinucleotides_profile$dinuc_Ref_rc,Dinucleotides_profile$dinuc_Alt,Dinucleotides_profile$dinuc_Alt_rc)

Dinucleotides_profile$dinuc_mutation_final <- paste0(Dinucleotides_profile$dinuc_Ref_final,">",Dinucleotides_profile$dinuc_Alt_final)

# All possible dinucleotide mutations (78)
a <- c("AA>CC","AA>CG","AA>CT","AA>GC","AA>GG","AA>GT","AA>TC","AA>TG","AA>TT",
       "AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
       "AG>CA","AG>CC","AG>CT","AG>GA","AG>GC","AG>GT","AG>TA","AG>TC","AG>TT",
       "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
       "CA>AC","CA>AG","CA>AT","CA>GC","CA>GG","CA>GT","CA>TC","CA>TG","CA>TT",
       "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
       "CG>AA","CG>AC","CG>AT","CG>GA","CG>GC","CG>TA",
       "GA>AC","GA>AG","GA>AT","GA>CC","GA>CG","GA>CT","GA>TC","GA>TG","GA>TT",
       "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
       "TA>AC","TA>AG","TA>AT","TA>CC","TA>CG","TA>GC")


Dinucleotides_template <- data.frame(a)
Dinucleotides_template$Ref <- substr(Dinucleotides_template[,1],1,2)
names(Dinucleotides_template) <- c("MutationType","Ref")
Dinucleotides_template10 <- data.frame("Ref"=c("AA","AC","AG","AT","CA","CC","CG","GA","GC","TA"),"Ref2"=c("AA","AC","AG","AT","CA","CC","CG","GA","GC","TA"))

# For dinucleotide number > 20 mutagens, Calculate cosimilarity
Dinucleotides_profile <- read.table("Dinucleotides.txt", sep = "\t", header = T, as.is = T)
a <- data.frame(table(Dinucleotides_profile$Sample.Name))
Dinucleotides_profile <- Dinucleotides_profile[Dinucleotides_profile$Sample.Name %in% a[a$Freq>=20,1],]

Dinucleotides_profile$firstmut <- ifelse(paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Alt)<=paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))),paste0(Dinucleotides_profile$Ref,">",Dinucleotides_profile$Alt),paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),">",as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))))
Dinucleotides_profile$secondmut <- ifelse(paste0(Dinucleotides_profile$Ref,Dinucleotides_profile$Alt)<=paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref_neighbor))),as.character(complement(DNAStringSet(Dinucleotides_profile$Alt_neighbor)))),paste0(Dinucleotides_profile$Ref_neighbor,">",Dinucleotides_profile$Alt_neighbor),paste0(as.character(complement(DNAStringSet(Dinucleotides_profile$Ref))),">",as.character(complement(DNAStringSet(Dinucleotides_profile$Alt)))))
Dinucleotides_profile$dimut <- paste0(Dinucleotides_profile$firstmut, ":",Dinucleotides_profile$secondmut)
sigfile_freq <- data.frame(table(Dinucleotides_profile$Sample.Name,Dinucleotides_profile$dimut))
names(sigfile_freq) <- c("Sample.Name","MutationType","Freq")
sigfile_freq <- merge(sigfile_freq,samples_details[,c("Sample.Name","Compound.Abbreviation","Concentration")],by="Sample.Name",all.x=T)
sigfile_freq$treatment <- paste0(sigfile_freq$Sample.Name, "_",sigfile_freq$Compound.Abbreviation,"_",sigfile_freq$Concentration)

control_sigset <-dcast(sigfile_freq,MutationType~treatment,value.var = "Freq")
control_sigset <- merge(dinuc_template,control_sigset,by="MutationType",all.x=T)
control_sigset[is.na(control_sigset)] <- 0
control_sigset_seletcted <- control_sigset

cossimil <- as.matrix(proxy::simil(t(control_sigset_seletcted[,4:dim(control_sigset_seletcted)[2]]), diag=FALSE, upper=FALSE, method="cosine",auto_convert_data_frames = FALSE))
cossimil[lower.tri(cossimil,diag=TRUE)]=NA
cossimil <- as.data.frame(as.table(cossimil))
cossimil=na.omit(cossimil)

write.table(cossimil,"Fig4D.txt",sep="\t",col.names = T, row.names = F, quote = F)
names(cossimil) <- c("sample1","sample2","simil")
cossimil$simil <- round(cossimil$simil,2)
filename=paste0("Fig4D.pdf")
pdf(file=filename, onefile=TRUE,width = 7.5,height = 6)
g1 <-ggplot(cossimil, aes(x=sample1, y=sample2)) + geom_tile(aes(fill=simil),colour="white")+geom_text(aes(label=paste(simil)),size=3)
g1 <-g1 +scale_fill_gradient2(midpoint=0.5,high="#B2182B", low="#2166AC",space="Lab",limits=c(0, 1))
g1 <-g1 +theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10),
               axis.text.y=element_text(size=10,colour = "black"),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
print(g1)
dev.off()





