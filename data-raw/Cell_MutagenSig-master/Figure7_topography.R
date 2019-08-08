dir.create("Figure7_results")
setwd("./Figure7_results/")
source("../Header.R")

replitime_length <- read.table("../00_data/MCF7_RepliTime_length.txt", sep = "\t", header = T, as.is = T)
transtrand_length <- read.table("../00_data/TranscriptionalStrand_length.txt", sep = "\t", header = T, as.is = T)
replitrand_length <- read.table("../00_data/ReplicativeStrand_length.txt", sep = "\t", header = T, as.is = T)

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

###################
# Fig7A  
# Transcriptional strand bias
###################
# 

sub_tab_all_info <- merge(sub_tab_all_info,samples_details,by="Sample.Name")


denovo_muts_child <- sub_tab_all_info
denovo_muts_child_bed <- denovo_muts_child[,c("Chrom","Pos","Pos","VariantID")]
denovo_muts_child_bed$Chrom <- paste0("chr",denovo_muts_child_bed$Chrom)
denovo_muts_child_bed[denovo_muts_child_bed$Chrom=="chr23","Chrom"] <- "chrX"
denovo_muts_child_bed[denovo_muts_child_bed$Chrom=="chr24","Chrom"] <- "chrY"
muts_control_mutagen_transstrand <- AddStrandInfo_intersect(sub_tab_all_info,"denovo_muts_control_mutagen_bed.txt","../00_data/TranscribStrand.uts.txt","../00_data/TranscribStrand.ts.txt","subs_control_mutagen","denovo_muts_control_mutagen_transstrand.txt")


# Remove background signature
# 53 treatments with stable signatures 
Sample_withSig <- read.table("../Figure3_results/Fig3B.txt",sep = "\t",header = T, as.is = T,quote = "\"")
sig_subs_full <- Sample_withSig[Sample_withSig$pvalue<=0.01 & Sample_withSig$profile_SNR>=2,]

muts_control_mutagen_transstrand <- read.table("./denovo_muts_control_mutagen_transstrand.txt",sep = "\t", header = T, as.is = T,quote = "\"")
  
muts_control_mutagen_transstrand2 <- muts_control_mutagen_transstrand[muts_control_mutagen_transstrand$Sample.Name %in% sig_subs_full$Sample.Name | muts_control_mutagen_transstrand$Group =="a_Control",]
muts_control_mutagen_transstrand2[muts_control_mutagen_transstrand2$Group=="a_Control",]$Sample <- paste0("Control","_",muts_control_mutagen_transstrand2[muts_control_mutagen_transstrand2$Group=="a_Control",]$Sample.Name)
denovo_6muttype_topography_byTreatment_control_mutagen(muts_control_mutagen_transstrand2,"Strand",transtrand_length,25,20,"FigS6_Transtrand_6muttype_control_mutagen")
muts_control_mutagen_transstrand3 <- muts_control_mutagen_transstrand[muts_control_mutagen_transstrand$Sample.Name %in% sig_subs_full$Sample.Name | muts_control_mutagen_transstrand$Group =="a_Control",]
muts_control_mutagen_transstrand3[muts_control_mutagen_transstrand3$Strand !="others","Strand"] <- paste0(muts_control_mutagen_transstrand3[muts_control_mutagen_transstrand3$Strand !="others",]$Ref2,">",muts_control_mutagen_transstrand3[muts_control_mutagen_transstrand3$Strand !="others",]$Alt2,"_",muts_control_mutagen_transstrand3[muts_control_mutagen_transstrand3$Strand !="others",]$Strand)
denovo_mutscount_sig_topography_strandbias(muts_control_mutagen_transstrand3,"Strand",transtrand_length,5,4)

# Chisq.test

  control_strand <- read.table("./FigS6_Transtrand_6muttype_control_mutagen.txt",sep = "\t",header = T, as.is = T)
  control_strand <- control_strand[control_strand$Sample.Name=="Control",c("Strand","Mutation","mean")]
  control_strand_dcast <- dcast(control_strand,Mutation~Strand,value.var="mean")
  control_strand_dcast$chisq.pvalue <- 1
  for(i in 1:6){
    control_strand_dcast[i,"chisq.pvalue"] <- chisq.test(c(control_strand_dcast[i,"-1"], control_strand_dcast[i,"1"]))$p.value
    
  }
  
  samplename_all <- NULL
  muttype_all <- NULL
  ratio_all <- NULL
  pvalue_all <- NULL
  for(i in 1:dim(sig_subs_full)[1]){
    mutagen_strand <- read.table(paste0("./Strand_",sig_subs_full[i,"Sample.Name"],"_",sig_subs_full[i,"Compound.Abbreviation"],"_exposure.txt"),sep = "\t",header = T, as.is = T)
    mutagen_strand <- mutagen_strand[,c("targetfeature","centroid","targetstrand")]
    names(mutagen_strand) <- c("targetfeature","mean","Strand")
    mutagen_strand$Mutation <- sub("\\_.*","",mutagen_strand$targetfeature)
    mutagen_strand$mean <- round(mutagen_strand$mean)
    mutagen_strand_dcast <- dcast(mutagen_strand[,c("mean","Strand","Mutation")],Mutation~Strand,value.var="mean")
    
    for(j in 1:dim(mutagen_strand_dcast)[1]){
      
      if(mutagen_strand_dcast[j,"-1"] >1 & mutagen_strand_dcast[j,"1"]>1 & sum(mutagen_strand_dcast[j,"-1"]+mutagen_strand_dcast[j,"1"])>=10){
        pvalue <- chisq.test(c(mutagen_strand_dcast[j,"-1"], mutagen_strand_dcast[j,"1"]))$p.value
        ratio <- mutagen_strand_dcast[j,2]/(mutagen_strand_dcast[j,2]+mutagen_strand_dcast[j,3])
      }else{
        pvalue <- 1
        ratio <- 0.5
      }
      
      pvalue_all <- c(pvalue_all,pvalue)
      ratio_all <- c(ratio_all,ratio)
      muttype_all <- c(muttype_all, mutagen_strand_dcast[j,1])
      samplename_all <- c(samplename_all,sig_subs_full[i,"Sample.Name"])
    }
  }
  chisq_result <- data.frame(samplename_all,muttype_all,ratio_all,pvalue_all)
  names(chisq_result) <- c("Sample.Name","Mutation","ratio","chisq_pvalue")
  chisq_result$P_adjust <- p.adjust(chisq_result$chisq_pvalue,method = "BH")
  chisq_result$flag <- ""
  chisq_result[chisq_result$P_adjust<=0.05,"flag"] <- "*"
  chisq_result[chisq_result$P_adjust<=0.01,"flag"] <- "**"
  chisq_result[chisq_result$P_adjust<=0.001,"flag"] <- "***"
  
  write.table(chisq_result,"transstrand_bias_chisq_adjust.txt",sep = "\t",col.names = T, row.names = F, quote = F)
  mutation_type_all_info <- merge(chisq_result,samples_details,by="Sample.Name")
  mutation_type_all_info$treatment <- paste0(mutation_type_all_info$Group,"_",mutation_type_all_info$Compound.Abbreviation,"_",mutation_type_all_info$Concentration,"_",mutation_type_all_info$Sample.Name)
  mutation_type_all_info <- mutation_type_all_info[order(mutation_type_all_info$Group,mutation_type_all_info$Compound.Abbreviation),]
  write.table(mutation_type_all_info,"Transtrand_6muttype_53treatment_control_chisq.txt",sep = "\t",col.names = T, row.names = F, quote = F)
  
  filename=paste0("Fig7A_Transtrand_6muttype_53treatment_control_chisq",".pdf")
  pdf(file=filename, onefile=TRUE,width = 8,height = 10)
  g1 <-ggplot(mutation_type_all_info, aes(x=Mutation, y=treatment)) + geom_tile(aes(fill=ratio),colour="black",size=0.25)
  g1 <- g1+geom_text(aes(label=paste(flag)),size=5)
  #  g1 <- g1+scale_y_discrete(limits=as.character(mutation_type_all_info$treatment))
  g1 <-g1 +scale_fill_gradient2(limits=c(0, 1),midpoint = 0.5)
  g1 <-g1 +theme(axis.text.x=element_text(size=10,colour = "black"),
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
  
  
  
##############################
#    Replicative strand bias
#    FigS5
##############################
AddRepliTimeInfo(sub_tab_all_info, "denovo_muts_control_mutagen_bed.txt", "../00_data/MCF7.all.compact.bed", "subs_ReplicatingTime.txt","muts_devono_replitime.txt")
muts_devono_replitime <- read.table("muts_devono_replitime.txt", sep = "\t", header = T, as.is = T, quote="\"")

Sample_withSig <- read.table("../Figure3_results/Fig3B.txt",sep = "\t",header = T, as.is = T,quote = "\"")
sig_subs_full <- Sample_withSig[Sample_withSig$pvalue<=0.01 & Sample_withSig$profile_SNR>=2,]
# Only do it for 53 signatures and control
denovo_mutscount_sig_topography2(muts_devono_replitime[(muts_devono_replitime$Sample.Name %in% sig_subs_full$Sample.Name | muts_devono_replitime$Group=="a_Control"),],"ReplicatingTime",replitime_length,8,4)



###########################################################
# Transcriptional strand bias in Replication timing regions
# For mutagens that have transcriptional strand bias, 
# compare the strand bias ratio in each replication timing regions
###########################################################
muts_devono_replitime <- read.table("muts_devono_replitime.txt", sep = "\t", header = T, as.is = T, quote="\"")
muts_devono_replitime2 <- muts_devono_replitime[muts_devono_replitime$ReplicatingTime!="others",]

RepliTime_list <- as.character(data.frame(table(muts_devono_replitime2[,"ReplicatingTime"]))[,1])
for(i in 1:length(RepliTime_list)){
  
  muts_curr <- muts_devono_replitime2[muts_devono_replitime2$ReplicatingTime==RepliTime_list[i],]
  muts_control_mutagen_transstrand <- AddStrandInfo_intersect(muts_curr,paste0(RepliTime_list[i],"_muts_curr_bed.txt"),"../00_data/TranscribStrand.uts.txt","../00_data/TranscribStrand.ts.txt",paste0(RepliTime_list[i],"_muts"),paste0(RepliTime_list[i],"_denovo_muts_control_mutagen_transstrand.txt"))
    
}
# Control
reptime_tranbias_all <- NULL
for(i in 1:length(RepliTime_list)){
  
  muts_curr <- muts_devono_replitime2[muts_devono_replitime2$ReplicatingTime==RepliTime_list[i],]
  muts_control_mutagen_transstrand <- read.table(paste0(RepliTime_list[i],"_denovo_muts_control_mutagen_transstrand.txt"),sep = "\t", header = T, as.is = T, quote="\"")
  reptime_tranbias <- denovo_6muttype_topography_byTreatment_control_aggragate(muts_control_mutagen_transstrand[muts_control_mutagen_transstrand$Group=="a_Control",],"Strand",transtrand_length,"NoPrint",3,4,paste0(RepliTime_list[i],"_","Transtrand_6muttype_control"))
  reptime_tranbias$ReplicatingTime <- RepliTime_list[i]
  reptime_tranbias_all <- rbind(reptime_tranbias_all,reptime_tranbias)
}
reptime_tranbias_all$ReplicatingTime <- as.numeric(reptime_tranbias_all$ReplicatingTime)
outputname=paste0("Control","_reptime_tranbias.pdf")
pdf(file=filename, onefile=TRUE,height=3,width = 12) 
d1 <- ggplot(reptime_tranbias_all,aes(x=ReplicatingTime,y=sum,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
d1 <- d1+ylab("mut. count")+theme(axis.title.x = element_text(size=15),
                                  axis.title.y = element_text(size=15),
                                  plot.title = element_text(size=10),
                                  axis.text.x=element_text(angle=90, vjust=0.5),
                                  panel.grid.minor.x=element_blank(),
                                  panel.grid.major.x=element_blank(),
                                  panel.grid.major.y = element_blank(),
                                  panel.grid.minor.y = element_blank(),
                                  panel.background = element_rect(fill = "white"),
                                  panel.border = element_rect(colour = "black", fill=NA))
d1 <- d1+facet_grid(.~ Mutation)
print(d1)
dev.off()

# Mutagen
muts_devono_replitime2 <- muts_devono_replitime2[muts_devono_replitime2$Group!="a_Control",]
mutagen_list <- c("MSM0.96","MSM0.13","MSM0.132","MSM0.12","MSM0.14","MSM0.103","MSM0.74","MSM0.42","MSM0.2","MSM0.92","MSM0.26",
                  "MSM0.55","MSM0.3","MSM0.27","MSM0.40","MSM0.93","MSM0.94","MSM0.95","MSM0.90","MSM0.41","MSM0.25","MSM0.83","MSM0.7","MSM0.84","MSM0.24")

for(j in 1:length(mutagen_list)){
  denovo_muts_child_rt_curr <- muts_devono_replitime2[muts_devono_replitime2$Sample.Name==mutagen_list[j],]
  reptime_tranbias_all <- NULL
  for(i in 1:length(RepliTime_list)){
    
    muts_control_mutagen_transstrand <- read.table(paste0(RepliTime_list[i],"_denovo_muts_control_mutagen_transstrand.txt"),sep = "\t", header = T, as.is = T, quote="\"")
    muts_control_mutagen_transstrand <- muts_control_mutagen_transstrand[muts_control_mutagen_transstrand$Sample.Name==mutagen_list[j],]
    reptime_tranbias <- denovo_6muttype_topography_byTreatment_control_aggragate(muts_control_mutagen_transstrand,"Strand",transtrand_length,"NoPrint",3,4,paste0(RepliTime_list[i],"_","Transtrand_6muttype_",muts_curr[1,"Treatment"],"_",mutagen_list[j]))
    reptime_tranbias$ReplicatingTime <- RepliTime_list[i]
    reptime_tranbias_all <- rbind(reptime_tranbias_all,reptime_tranbias)
  }
  reptime_tranbias_all$ReplicatingTime <- as.numeric(reptime_tranbias_all$ReplicatingTime)
  outputname=paste0(mutagen_list[j],"_",denovo_muts_child_rt_curr[1,"Treatment"],"_reptime_tranbias.pdf")
  pdf(file=outputname, onefile=TRUE,height=3,width = 12) 
  d1 <- ggplot(reptime_tranbias_all,aes(x=ReplicatingTime,y=sum,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
  d1 <- d1+ylab("mut. count")+theme(axis.title.x = element_text(size=15),
                                    axis.title.y = element_text(size=15),
                                    plot.title = element_text(size=10),
                                    axis.text.x=element_text(angle=90, vjust=0.5),
                                    panel.grid.minor.x=element_blank(),
                                    panel.grid.major.x=element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    panel.background = element_rect(fill = "white"),
                                    panel.border = element_rect(colour = "black", fill=NA))
  d1 <- d1+facet_grid(.~ Mutation)
  print(d1)
  dev.off()
  
}






  if(TRUE){
    Tab2Bed(sub_tab_all_info,"denovo_muts_bed.txt")
    
    # nohup intersectBed -a denovo_muts_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/MCF7.all.compact.bed -wo > subs_ReplicatingTime.txt &
    sub_tab_all_info <- merge(sub_tab_all_info,samples_details,by="Sample.Name")
    
    denovo_muts_child <- sub_tab_all_info
    subs_ReplicatingTime <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/13_subs_newnames/StrandBias/ReplicatingTime/subs_ReplicatingTime.txt",sep = "\t",header = F,as.is = T)
    subs_ReplicatingTime_short <- subs_ReplicatingTime[,c(4,9)]
    names(subs_ReplicatingTime_short) <- c("VariantID","ReplicatingTime")
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==1,"ReplicatingTime"] <- "j"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==2,"ReplicatingTime"] <- "i"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==3,"ReplicatingTime"] <- "h"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==4,"ReplicatingTime"] <- "g"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==5,"ReplicatingTime"] <- "f"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==6,"ReplicatingTime"] <- "e"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==7,"ReplicatingTime"] <- "d"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==8,"ReplicatingTime"] <- "c"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==9,"ReplicatingTime"] <- "b"
    subs_ReplicatingTime_short[subs_ReplicatingTime_short$ReplicatingTime==10,"ReplicatingTime"] <- "a"
    
    denovo_muts_child_rt <- merge(denovo_muts_child,subs_ReplicatingTime_short,by="VariantID",all.x=T)
    denovo_muts_child_rt[is.na(denovo_muts_child_rt)] <- "others"
    denovo_muts_child_rt <- denovo_muts_child_rt[denovo_muts_child_rt$ReplicatingTime!="others",]
    RepliTime_list <- as.character(data.frame(table(denovo_muts_child_rt[,"ReplicatingTime"]))[,1])
    
    for(i in 1:length(RepliTime_list)){
      
      muts_curr <- denovo_muts_child_rt[denovo_muts_child_rt$ReplicatingTime==RepliTime_list[i],]
      
      # control + mutagen
      Tab2Bed(muts_curr,paste0(RepliTime_list[i],"_muts_curr_bed.txt"))
    }
    
    # nohup intersectBed -a a_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > a_muts_curr_uts.txt &
    # nohup intersectBed -a a_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > a_muts_curr_ts.txt &
    #  nohup intersectBed -a b_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > b_muts_curr_uts.txt &
    #  nohup intersectBed -a b_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > b_muts_curr_ts.txt &
    #  nohup intersectBed -a c_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > c_muts_curr_uts.txt &
    #  nohup intersectBed -a c_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > c_muts_curr_ts.txt &
    #  nohup intersectBed -a d_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > d_muts_curr_uts.txt &
    #  nohup intersectBed -a d_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > d_muts_curr_ts.txt &
    #  nohup intersectBed -a e_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > e_muts_curr_uts.txt &
    #  nohup intersectBed -a e_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > e_muts_curr_ts.txt &
    #  nohup intersectBed -a f_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > f_muts_curr_uts.txt &
    #  nohup intersectBed -a f_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > f_muts_curr_ts.txt &
    #  nohup intersectBed -a g_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > g_muts_curr_uts.txt &
    #  nohup intersectBed -a g_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > g_muts_curr_ts.txt &
    #  nohup intersectBed -a h_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > h_muts_curr_uts.txt &
    #  nohup intersectBed -a h_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > h_muts_curr_ts.txt &
    #  nohup intersectBed -a i_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > i_muts_curr_uts.txt &
    #  nohup intersectBed -a i_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > i_muts_curr_ts.txt &
    #  nohup intersectBed -a j_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.uts.txt -wo > j_muts_curr_uts.txt &
    #  nohup intersectBed -a j_muts_curr_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/TranscribStrand.ts.txt -wo > j_muts_curr_ts.txt &
    
    
    # Control
    reptime_tranbias_all <- NULL
    for(i in 1:length(RepliTime_list)){
      
      muts_curr <- denovo_muts_child_rt[denovo_muts_child_rt$ReplicatingTime==RepliTime_list[i],]
      muts_control_mutagen_transstrand <- AddStrandInfo(paste0(RepliTime_list[i],"_muts_curr_uts.txt"),paste0(RepliTime_list[i],"_muts_curr_ts.txt"), muts_curr, paste0(RepliTime_list[i],"_muts_curr_transstrand.txt"))
      reptime_tranbias <- denovo_6muttype_topography_byTreatment_control_aggragate(muts_control_mutagen_transstrand[muts_control_mutagen_transstrand$Group=="a_Control",],"Strand",transtrand_length,"NoPrint",3,4,paste0(RepliTime_list[i],"_","Transtrand_6muttype_control"))
      reptime_tranbias$ReplicatingTime <- RepliTime_list[i]
      reptime_tranbias_all <- rbind(reptime_tranbias_all,reptime_tranbias)
    }
    outputname=paste0("Control","_reptime_tranbias.pdf")
    filename=paste0(outputname,".pdf")
    pdf(file=filename, onefile=TRUE,height=3,width = 12) 
    d1 <- ggplot(reptime_tranbias_all,aes(x=ReplicatingTime,y=sum,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
    d1 <- d1+ylab("mut. count")+theme(axis.title.x = element_text(size=15),
                                      axis.title.y = element_text(size=15),
                                      plot.title = element_text(size=10),
                                      axis.text.x=element_text(angle=90, vjust=0.5),
                                      panel.grid.minor.x=element_blank(),
                                      panel.grid.major.x=element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.minor.y = element_blank(),
                                      panel.background = element_rect(fill = "white"),
                                      panel.border = element_rect(colour = "black", fill=NA))
    d1 <- d1+facet_grid(.~ Mutation)
    print(d1)
    dev.off()
    
    
    # Mutagen
    denovo_muts_child_rt <- denovo_muts_child_rt[denovo_muts_child_rt$Group!="a_Control",]
    #mutagen_list <- data.frame(table(denovo_muts_child_rt$Sample.Name))
    mutagen_list <- c("MSM0.96","MSM0.13","MSM0.132","MSM0.12","MSM0.14","MSM0.103","MSM0.74","MSM0.42","MSM0.2","MSM0.92","MSM0.92","MSM0.26",
                      "MSM0.55","MSM0.3","MSM0.27","MSM0.40","MSM0.93","MSM0.94","MSM0.95","MSM0.90","MSM0.41","MSM0.25","MSM0.83","MSM0.7","MSM0.84","MSM0.24")
    
    for(j in 1:length(mutagen_list)){
      denovo_muts_child_rt_curr <- denovo_muts_child_rt[denovo_muts_child_rt$Sample.Name==mutagen_list[j],]
      reptime_tranbias_all <- NULL
      for(i in 1:length(RepliTime_list)){
        
        muts_curr <- denovo_muts_child_rt_curr[denovo_muts_child_rt_curr$ReplicatingTime==RepliTime_list[i],]
        muts_control_mutagen_transstrand <- AddStrandInfo(paste0(RepliTime_list[i],"_muts_curr_uts.txt"),paste0(RepliTime_list[i],"_muts_curr_ts.txt"), muts_curr, paste0(RepliTime_list[i],"_muts_curr_transstrand.txt"))
        reptime_tranbias <- denovo_6muttype_topography_byTreatment_control_aggragate(muts_control_mutagen_transstrand,"Strand",transtrand_length,"NoPrint",3,4,paste0(RepliTime_list[i],"_","Transtrand_6muttype_",muts_curr[1,"Treatment"],"_",mutagen_list[j]))
        reptime_tranbias$ReplicatingTime <- RepliTime_list[i]
        reptime_tranbias_all <- rbind(reptime_tranbias_all,reptime_tranbias)
      }
      outputname=paste0(mutagen_list[j],"_",denovo_muts_child_rt_curr[1,"Treatment"],"_reptime_tranbias.pdf")
      filename=paste0(outputname,".pdf")
      pdf(file=filename, onefile=TRUE,height=3,width = 12) 
      d1 <- ggplot(reptime_tranbias_all,aes(x=ReplicatingTime,y=sum,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
      d1 <- d1+ylab("mut. count")+theme(axis.title.x = element_text(size=15),
                                        axis.title.y = element_text(size=15),
                                        plot.title = element_text(size=10),
                                        axis.text.x=element_text(angle=90, vjust=0.5),
                                        panel.grid.minor.x=element_blank(),
                                        panel.grid.major.x=element_blank(),
                                        panel.grid.major.y = element_blank(),
                                        panel.grid.minor.y = element_blank(),
                                        panel.background = element_rect(fill = "white"),
                                        panel.border = element_rect(colour = "black", fill=NA))
      d1 <- d1+facet_grid(.~ Mutation)
      print(d1)
      dev.off()
      
      
    }
    
  }
  


