library(Biostrings)
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(reshape2)
library(proxy)
library(stats)
library(gplots)
library(RColorBrewer)
library(dndscv)
library("VariantAnnotation")
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)


gen_muttype_new <- function(CTsubs){
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/14_subs_0913/MutationType_template.txt", sep = "\t", header = T, as.is = T)
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs$Sample,CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}
plotCountbasis_average_se <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,2:dim(muts_basis)[2]])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("MutationType"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N),snr=abs(mean/sd))
  muts_basis_melt_summary$mutation <- substr(muts_basis_melt_summary$MutationType,3,5)
  muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$mutation),]
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=mean,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),size=.2,width=0.5) #+scale_y_continuous(limits=c(0, 350),breaks=seq(0, 350, 100))
  p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4,colour = "black"),
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
  #  return(muts_basis_melt_summary)
  write.table(muts_basis_melt_summary,paste0(outputname, ".txt"), sep = "\t", col.names = T, row.names = F, quote = F)
}

sd_highD <- function(target_matrix, direction){
  
  centroid <- rowMeans(target_matrix)
  target_matrix_diff <- target_matrix-centroid
  
  #### Project target_matrix_diff to given direction
  # project_vector <- apply(target_matrix_diff,2,function(x) project(x, direction, type='length'))
  cos_vector_square <- apply(target_matrix_diff,2,function(x) (norm(as.matrix(x),"f")*cos_similarity(x,direction))^2)
  sd_matrix <- sqrt(sum(cos_vector_square)/(dim(target_matrix_diff)[2]-1))
  #sd_matrix <- 
  return(sd_matrix)
  
}
bootstrapGenomesfun <- function(genomes){
  
  return(apply(genomes, 2, function(x) rmultinom(1, sum(x), x)))
}
cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
}

ExtractSig_centroid_subs_noerrorbar <- function(control_profile_set, compound_profile_set, sampling_number,boundary,muttype_template,h,w,outputname) {
  
  # Test convergence of sig
  if(TRUE){
    result_profile <- NULL
    for(j in 1:dim(compound_profile_set)[2]){
      a <- RemoveBackground_single(control_profile_set,compound_profile_set[,j],1.65)
      result_profile <- cbind(result_profile,a)
    }
    
    min_b <- 1
    max_b <- 0
    for(j in 1:dim(compound_profile_set)[2]){
      if((j+1)<=dim(compound_profile_set)[2]){
        for(k in (j+1):dim(compound_profile_set)[2]){
          b <- cos_similarity(result_profile[,j],result_profile[,k])
          min_b <- ifelse(b<min_b,b,min_b)
          max_b <- ifelse(b>max_b,b,max_b)
          print(paste0("min_b:",min_b,"; max_b:",max_b))
        }
      }
    }
    stability_sig <- c(min_b,max_b)
    
  }
  
  
  compound_sig <- RemoveBackground_centroid(control_profile_set,compound_profile_set,1000,1.65)
  compound_sig$MutationType <- rownames(compound_sig)
  muts_basis_melt_summary <- merge(muttype_template, compound_sig,by="MutationType",all.x=T)
  muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
  muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$Mutation),]
  
  muts_basis_melt_summary$percentage <- muts_basis_melt_summary[,"centroid"]/sum(muts_basis_melt_summary[,"centroid"])
  muts_basis_melt_summary$percentage_sd <- muts_basis_melt_summary[,"sd"]/sum(muts_basis_melt_summary[,"centroid"])
  
  write.table(muts_basis_melt_summary,paste0(outputname, "_exposure.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  filename <- paste0(outputname, "_Signature.pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=centroid,fill=Mutation))+ geom_bar(stat="identity",position="dodge", width=.7)+xlab("Substitution Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(paste0(outputname,"_exposure"))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4.5,colour = "black"),
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
  
  
  o <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=percentage,fill=Mutation))+ geom_bar(stat="identity",position="dodge", width=.7)+xlab("Substitution Types")+ylab("Percentage")
  o <- o+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(paste0(outputname,"_signature","(stability=",round(max_b,2),")"))
  o <- o+scale_y_continuous(labels = scales::percent,limits=c(0,0.25),breaks=(seq(0,0.25,0.1)))
  o <- o+scale_fill_manual(values=mypalette)
  o <- o+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4.5,colour = "black"),
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
  
  grid.arrange(p,o,ncol=1)
  dev.off()
  
  return(stability_sig)
  
}
RemoveBackground_single <- function(background_profile_set, sig_profile,boundary){
  
  # Range of background
  centroid_background <- rowMeans(background_profile_set)
  sd_background <- apply(background_profile_set,1,sd)
  boundary_background <- centroid_background+boundary*sd_background
  
  
  diff_all_boundary <- sig_profile-boundary_background
  diff_all <- sig_profile-centroid_background
  
  diff_all[which(diff_all_boundary<0)] <- 0
  return(diff_all)
  
}
RemoveBackground_centroid <- function(background_profile_set, sig_profile_set, sampling_number,boundary){
  
  # Range of background
  centroid_background <- rowMeans(background_profile_set)
  sd_background <- apply(background_profile_set,1,sd)
  boundary_background <- centroid_background+boundary*sd_background
  
  # Bootstrap subclones 
  diff_mean_all <- NULL
  for(bt_num in 1:sampling_number){
    
    RepCompound <- sig_profile_set[,sample(1:dim(sig_profile_set)[2],dim(sig_profile_set)[2],replace = T)]
    bootstrapCompound <- bootstrapGenomesfun(RepCompound)
    
    diff_all_boundary <- rowMeans(bootstrapCompound)-boundary_background
    diff_all <- rowMeans(bootstrapCompound)-centroid_background
    
    diff_all[which(diff_all_boundary<0)] <- 0
    diff_mean_all <- cbind(diff_mean_all,diff_all)
    
  }
  
  diff_centroid_all <- rowMeans(diff_mean_all)
  diff_quantile_all_sd <- apply(diff_mean_all,1,sd)
  
  pure_sig <- data.frame(cbind(diff_centroid_all,diff_quantile_all_sd))
  names(pure_sig) <- c("centroid","sd")
  
  return(pure_sig)
  
}


Mutagenicity_ratio <- function(background_profile_set, sig_profile_set){
  
  # Bootstrap subclones 
  MutagenicRatio_all <- NULL
  diff_mean_all <- NULL
  for(i in 1:length(background_profile_set)){
    for(j in 1:length(sig_profile_set)){
      MutagenicRatio <- (sig_profile_set[j]-background_profile_set[i])/background_profile_set[i]
      MutagenicRatio_all <- c(MutagenicRatio_all,MutagenicRatio)
      
    }
    
  }
  
  return(c(mean(MutagenicRatio_all),sd(MutagenicRatio_all)/sqrt(length(MutagenicRatio_all))))
  
}



#####################################
#
# Figure 4 Double substitutions
#
#####################################
gen_muttype_Dinucleotides <- function(CTsubs,SelectCol){
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/12_subs/Dinucleiotides/Dinucleotides_template.txt", sep = "\t", header = T, as.is = T)
  
  sigfile_freq <- data.frame(table(CTsubs[,SelectCol],CTsubs$dinuc_mutation_final))
  names(sigfile_freq) <- c("SelectPara","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~SelectPara,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  return(control_sigset)
  
}
plot_allsample_Dinucleotides_profile <- function(muttype_catalogue,colnum,h,w,outputname){
  
  muts_basis_melt <- melt(muttype_catalogue,id=c("MutationType","Ref"))
  names(muts_basis_melt) <- c("MutationType","Ref","sample","count")
  
  mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe")
  
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=Ref,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("Dinucleotide mutation type")+ylab("")
  #  p <- p+coord_cartesian(ylim=c(0, max(muttype_freq[,"freq"])))
  p <- p+scale_x_discrete(limits = as.character(muttype_catalogue$MutationType))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black"),
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
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  print(p)
  dev.off()
} 
FindDinucleotides <- function(allsubs) {
  dinuclist <- NULL
  samplelist <- data.frame(table(allsubs$Sample))
  for(i in 1:dim(samplelist)[1]){
    print(i)
    samplesubs <- allsubs[allsubs$Sample==samplelist[i,1],]
    a <- samplesubs$Pos
    b <- samplesubs$Ref
    c <- samplesubs$Alt
    samplesubs$Pos_neighbor <- c(a[-1],0)
    samplesubs$neigbor_dist <- samplesubs$Pos-samplesubs$Pos_neighbor
    samplesubs$Ref_neighbor <- c(b[-1],"N")
    samplesubs$Alt_neighbor <- c(c[-1],"N")
    
    dinuc_index <- which(samplesubs$neigbor_dist==-1)
    print(length(dinuc_index))
    if(length(dinuc_index)>0){
      dinuclist <- rbind(dinuclist,samplesubs[dinuc_index,])
      dinuclist <- rbind(dinuclist,samplesubs[dinuc_index+1,])
    }
  }
  return(dinuclist)
}

#######################
#
# Figure 5 Indels
#
#######################
plotbasis_aggragated_indel_29types_6 <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_basis$aggragate <- rowSums(muts_basis[,1:(dim(muts_basis)[2]-1)])/sum(muts_basis[,1:(dim(muts_basis)[2]-1)])
  
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/29_indel_29channels_newnames/indel_29channels_6.txt",sep = "\t",header = T, as.is = T)
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  
  # E69F00(light orange, +C); D55E00(dark orange, +T); F0E442(yellow, +>1); 56B4E9(light blue, -C); 0072B2(dark blue, -T); 009E73(green, ->1); CC79A7(pink, MhDel)
  #indel_mypalette <- c("#009E73","#56B4E9","#0072B2","#F0E442","#E69F00","#D55E00")
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  indel_mypalette <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00","#911eb4")
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LongRep",
                       "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LongRep",
                       "[+>1]NonRep", "[+>1]Rep",
                       "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LongRep",
                       "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LongRep",
                       "[->1]Others","[->1]Rep","[->1]Mh") 
  
  indel_labels <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LR",
                    "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LR",
                    "[+>1]NonR", "[+>1]Rep",
                    "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LR",
                    "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LR",
                    "[->1]Oth.","[->1]Rep","Mh") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=aggragate,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  p <- p+scale_y_continuous(limits=c(0,0.4),breaks=(seq(0,0.4,0.1)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
}
plotCountbasis_aggragated_indel_29types_6 <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_basis$aggragate <- rowSums(muts_basis[,1:(dim(muts_basis)[2]-1)])
  
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/29_indel_29channels_newnames/indel_29channels_6.txt",sep = "\t",header = T, as.is = T)
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  
  # E69F00(light orange, +C); D55E00(dark orange, +T); F0E442(yellow, +>1); 56B4E9(light blue, -C); 0072B2(dark blue, -T); 009E73(green, ->1); CC79A7(pink, MhDel)
  #indel_mypalette <- c("#009E73","#56B4E9","#0072B2","#F0E442","#E69F00","#D55E00")
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  indel_mypalette <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00","#911eb4")
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LongRep",
                       "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LongRep",
                       "[+>1]NonRep", "[+>1]Rep",
                       "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LongRep",
                       "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LongRep",
                       "[->1]Others","[->1]Rep","[->1]Mh") 
  
  indel_labels <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LR",
                    "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LR",
                    "[+>1]NonR", "[+>1]Rep",
                    "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LR",
                    "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LR",
                    "[->1]Oth.","[->1]Rep","Mh") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=aggragate,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}
ExtractSig_aggregate_indels <- function(control_profile_set, compound_profile_set, boundary,muttype_template,h,w,outputname) {
  
  # Test convergence of sig
  if(TRUE){
    result_profile <- NULL
    for(j in 1:dim(compound_profile_set)[2]){
      a <- RemoveBackground_single(control_profile_set,compound_profile_set[,j],boundary)
      result_profile <- cbind(result_profile,a)
    }
    
    min_b <- 1
    max_b <- 0
    for(j in 1:dim(compound_profile_set)[2]){
      if((j+1)<=dim(compound_profile_set)[2]){
        for(k in (j+1):dim(compound_profile_set)[2]){
          b <- cos_similarity(result_profile[,j],result_profile[,k])
          min_b <- ifelse(b<min_b,b,min_b)
          max_b <- ifelse(b>max_b,b,max_b)
          print(paste0("min_b:",min_b,"; max_b:",max_b))
        }
      }
    }
    stability_sig <- c(min_b,max_b)
    
  }
  
  
  compound_aggregate <- rowMeans(compound_profile_set)
  compound_sig <- as.data.frame(RemoveBackground_single(control_profile_set,compound_aggregate,boundary))
  names(compound_sig) <- "exposure"
  compound_sig$indelsubtype <- rownames(compound_sig)
  muts_basis_melt_summary <- merge(muttype_template, compound_sig,by="indelsubtype",all.x=T)
  muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
  muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$indeltype),]
  
  muts_basis_melt_summary$percentage <- muts_basis_melt_summary[,"exposure"]/sum(muts_basis_melt_summary[,"exposure"])
  
  write.table(muts_basis_melt_summary,paste0(outputname, "_exposure.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  indel_mypalette <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00","#911eb4")
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LongRep",
                       "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LongRep",
                       "[+>1]NonRep", "[+>1]Rep",
                       "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LongRep",
                       "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LongRep",
                       "[->1]Others","[->1]Rep","[->1]Mh") 
  
  indel_labels <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LR",
                    "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LR",
                    "[+>1]NonR", "[+>1]Rep",
                    "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LR",
                    "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LR",
                    "[->1]Oth.","[->1]Rep","[-]Mh") 
  
  
  filename <- paste0(outputname, "_IndelSignature.pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=indelsubtype, y=exposure,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = indel_positions, labels =  indel_labels)+ggtitle(paste0(outputname,"_exposure"))
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black"),
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
  
  
  o <- ggplot(data=muts_basis_melt_summary, aes(x=indelsubtype, y=percentage,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  o <- o+scale_y_continuous(labels = scales::percent) #limits=c(0,0.3),breaks=(seq(0,1,0.3)),
  o <- o+scale_x_discrete(limits = indel_positions, labels =  indel_labels)+ggtitle(paste0(outputname,"_signature"))
  o <- o+scale_fill_manual(values=indel_mypalette)
  o <- o+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black"),
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
  
  grid.arrange(p,o,ncol=1)
  dev.off()
  
  return(stability_sig)
}

#######################
#
# Figure 7 Topography
#
#######################
# bed file is 0-based, half-closed-half-open 
Tab2Bed <- function(muts, outputname){ # bed file is 0-based, half-closed-half-open 
  muts_bed <- muts[,c("Chrom","Pos","Pos","VariantID")]
  muts_bed$Chrom <- paste0("chr",muts_bed$Chrom)
  muts_bed[muts_bed$Chrom=="chr23","Chrom"] <- "chrX"
  muts_bed[muts_bed$Chrom=="chr24","Chrom"] <- "chrY"
  muts_bed[,2] <- muts_bed[,2]-1
  write.table(muts_bed,outputname,sep="\t",col.names = F, row.names = F, quote = F)
}
AddStrandInfo <- function(mutfile1, mutfile2,muts_context,outputname){
  mut_strand1 <- read.table(mutfile1,sep = "\t",header = F,as.is = T)
  mut_strand2 <- read.table(mutfile2,sep = "\t",header = F,as.is = T)
  mut_strand <- rbind(mut_strand1,mut_strand2)
  mut_strand_short <- mut_strand[,c(4,8)]
  names(mut_strand_short) <- c("VariantID","Strand_original")
  
  muts_withstrandinfo <- merge(muts_context,mut_strand_short,by="VariantID",all.x=T)
  muts_withstrandinfo[is.na(muts_withstrandinfo)] <- "others"
  muts_withstrandinfo[(muts_withstrandinfo$Strand_original == "Leading" | muts_withstrandinfo$Strand_original == 1),]$Strand_original <- 1
  muts_withstrandinfo[(muts_withstrandinfo$Strand_original == "Lagging" | muts_withstrandinfo$Strand_original == -1),]$Strand_original <- -1
  
  muts_withstrandinfo$Strand <- muts_withstrandinfo$Strand_original
  
  CTsubs_copy <- muts_withstrandinfo
  CTsubs <- muts_withstrandinfo
  CTsubs$Alt2 <- CTsubs$Alt
  CTsubs$Ref2 <- CTsubs$Ref
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt2 <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref2 <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs[(CTsubs$Ref2!=CTsubs$Ref & CTsubs$Strand_original=="1"),]$Strand <- -1
  CTsubs[(CTsubs$Ref2!=CTsubs$Ref & CTsubs$Strand_original=="-1"),]$Strand <- 1
  
  
  write.table(CTsubs,outputname,sep="\t",col.names = T, row.names = F, quote = F)
  return(CTsubs)
}

AddRepliTimeInfo <- function(mutlist, mutBedfile,featureBedfile,intersectResultfile,outputfilename){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile," -wo > ", intersectResultfile)
  
  #"nohup intersectBed -a denovo_muts_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/MCF7.all.compact.bed -wo > subs_ReplicatingTime.txt &"
  try(system(intersectBed_command))
  
  subs_ReplicatingTime <- read.table(intersectResultfile,sep = "\t",header = F,as.is = T)
  subs_ReplicatingTime_short <- subs_ReplicatingTime[,c(4,9)]
  names(subs_ReplicatingTime_short) <- c("VariantID","ReplicatingTime")
  subs_ReplicatingTime_short$ReplicatingTime=abs(subs_ReplicatingTime_short$ReplicatingTime-11)
  denovo_muts_rt <- merge(mutlist,subs_ReplicatingTime_short,by="VariantID",all.x=T)
  denovo_muts_rt[is.na(denovo_muts_rt)] <- "others"
  write.table(denovo_muts_rt,outputfilename,sep = "\t", col.names = T, row.names = F, quote = F)
  
}
AddStrandInfo_intersect <- function(mutlist, mutBedfile,featureBedfile_strand1,featureBedfile_strand2,intersectResultfile,outputfilename){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command_1 <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile_strand1," -wo > ", paste0(intersectResultfile,"_1.txt"))
  intersectBed_command_2 <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile_strand2," -wo > ", paste0(intersectResultfile,"_2.txt"))
  
  # nohup intersectBed -a denovo_muts_control_mutagen_bed.txt -b ./00_data/TranscribStrand.uts.txt -wo > subs_control_mutagen_uts.txt &
  # nohup intersectBed -a denovo_muts_control_mutagen_bed.txt -b ./00_data/TranscribStrand.ts.txt -wo > subs_control_mutagen_ts.txt &
  try(system(intersectBed_command_1))
  try(system(intersectBed_command_2))
  
  AddStrandInfo(paste0(intersectResultfile,"_1.txt"),paste0(intersectResultfile,"_2.txt"),mutlist,outputfilename)
}
denovo_6muttype_topography_byTreatment_control_mutagen <- function(denovo_muts_child_feature, feature,featurelength,h,w, outputname){
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  featurelength$all <-rowSums(featurelength[,2:5])
  
  denovo_muts_child_feature <- denovo_muts_child_feature[denovo_muts_child_feature[,feature] != "others",]
  denovo_muts_child_feature$Mutation <- paste0(denovo_muts_child_feature$Ref2,">",denovo_muts_child_feature$Alt2)
  denovo_muts_child_dis <- data.frame(table(denovo_muts_child_feature[,"Sample"],denovo_muts_child_feature[,feature], denovo_muts_child_feature[,"Mutation"]))
  names(denovo_muts_child_dis) <- c("Sample",feature,"Mutation", "Freq")
  
  denovo_muts_child_dis$Sample.Name <- sub("\\_.*","",denovo_muts_child_dis$Sample)
  denovo_muts_child_dis[denovo_muts_child_dis$Group=="Control","Sample.Name"]="Control"
  
  gtc <- ddply(denovo_muts_child_dis, c(feature,"Mutation","Sample.Name"),summarise, N=length(Sample.Name),mean=mean(Freq),sd=sd(Freq),se=sd/sqrt(N))
  gtc$targetfeature <- gtc[,feature]
  gtc_density <- merge(gtc,featurelength,by=feature,all.x=T)
  gtc_density$mean_d <- gtc_density$mean/gtc_density$all*1000000
  gtc_density$se_d <- gtc_density$se/gtc_density$all*1000000
  gtc_density <- merge(gtc_density,samples_details,by="Sample.Name",all.x=T)
  gtc_density$treatment <- paste0(gtc_density$Sample.Name,"_",gtc_density$Compound.Abbreviation,"_",gtc_density$Concentration)
  featurelist <- data.frame(table(gtc[,feature]))[,1]
  filename=paste0(outputname,".pdf")
  pdf(file=filename, onefile=TRUE,height=h,width = w) 
  d1 <- ggplot(gtc_density,aes(x=Mutation,y=mean_d,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
  d1 <- d1+geom_errorbar(aes(ymin=mean_d-se_d,ymax=mean_d+se_d),position=position_dodge(.9),width=.2)
  d1 <- d1+ylab("mut. density")+theme(axis.title.x = element_text(size=15),
                                      axis.title.y = element_text(size=15),
                                      plot.title = element_text(size=10),
                                      axis.text.x=element_text(angle=90, vjust=0.5),
                                      panel.grid.minor.x=element_blank(),
                                      panel.grid.major.x=element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.minor.y = element_blank(),
                                      panel.background = element_rect(fill = "white"),
                                      panel.border = element_rect(colour = "black", fill=NA))
  d1 <- d1+facet_wrap(~treatment,ncol=7,scales="free")
  print(d1)
  dev.off()
  
  write.table(gtc_density,paste0(outputname,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
}

# Use ExactSig function, count for Strand Bias
denovo_mutscount_sig_topography_strandbias <- function(denovo_muts_child_feature, feature,featurelength,h,w){
  featurelength$all <-rowSums(featurelength[,2:5])
  
  denovo_muts_child_feature <- denovo_muts_child_feature[denovo_muts_child_feature[,feature] != "others",]
  denovo_muts_child_dis <- data.frame(table(denovo_muts_child_feature[,"Sample"],denovo_muts_child_feature[,feature]))
  names(denovo_muts_child_dis) <- c("Sample",feature, "Freq")
  denovo_muts_child_dis$Sample.Name <- sub("\\_.*","",denovo_muts_child_dis$Sample)
  denovo_muts_child_dis <- merge(denovo_muts_child_dis,samples_details,by="Sample.Name")
  denovo_muts_child_dis$Treatment <- paste0(denovo_muts_child_dis$Sample.Name,"_",denovo_muts_child_dis$Compound.Abbreviation,"_",denovo_muts_child_dis$Concentration)
  denovo_muts_child_dis$targetfeature <- denovo_muts_child_dis[,feature]
  # Control
  control_muts <- denovo_muts_child_dis[denovo_muts_child_dis$Group=="a_Control",]
  control_profile <- dcast(control_muts[,c("Sample","targetfeature","Freq")], targetfeature~Sample, value.var="Freq")
  rownames(control_profile) <- control_profile$targetfeature
  
  
  # Compounds
  allcompound_profile <- denovo_muts_child_dis[denovo_muts_child_dis$Group!="a_Control",]
  compoundlist <- data.frame(table(allcompound_profile$Sample.Name))
  names(compoundlist) <- c("Sample.Name","Freq")
  #compoundlist$P_occurance <- 0
  mut_template <- data.frame(control_profile$targetfeature)
  names(mut_template) <- "targetfeature"
  mutation_type_all <- NULL
  for(i in 1:dim(compoundlist)[1]){
    print(i)
    currentcompound <- denovo_muts_child_dis[denovo_muts_child_dis$Sample.Name==compoundlist[i,"Sample.Name"],]
    compound_profile <- dcast(currentcompound[,c("Sample","targetfeature","Freq")], targetfeature~Sample, value.var="Freq")
    rownames(compound_profile) <- compound_profile$targetfeature
    #  ExtractSig(control_profile[,-1],compound_profile[,-1],1000,1.65,mut_template,feature,4.2,7,paste0(feature,as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Compound.Abbreviation"]))
    compound_sig <- RemoveBackground_centroid(control_profile[,-1],compound_profile[,-1],1000,1.65)
    compound_sig$targetfeature <- rownames(compound_sig)
    muts_basis_melt_summary <- merge(mut_template, compound_sig,by="targetfeature",all.x=T)
    muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
    muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$targetfeature),]
    
    muts_basis_melt_summary$percentage <- muts_basis_melt_summary[,"centroid"]/sum(muts_basis_melt_summary[,"centroid"])
    muts_basis_melt_summary$percentage_sd <- muts_basis_melt_summary[,"sd"]/sum(muts_basis_melt_summary[,"centroid"])
    muts_basis_melt_summary$targetstrand <- substr(muts_basis_melt_summary$targetfeature,5,6)
    outputname=paste0(feature,"_",as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Compound.Abbreviation"])
    
    write.table(muts_basis_melt_summary,paste0(outputname, "_exposure.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
    
    #mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#808080")
    if(FALSE){
      filename <- paste0(outputname, "_Signature.pdf")
      pdf(file=filename, onefile=TRUE,width=w,height=h)
      p <- ggplot(data=muts_basis_melt_summary, aes(x=targetfeature, y=centroid,fill=targetstrand))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab(feature)+ylab("Count")
      p <- p+geom_errorbar(aes(ymin=centroid-sd,ymax=centroid+sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)#+scale_y_continuous(limits=c(-1,1),breaks=(seq(-1,1,0.2)),labels = scales::percent)
      p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$targetfeature))+ggtitle(paste0(outputname,"_exposure"))
      #   p <- p+scale_fill_manual(values=mypalette)
      p <- p+theme(axis.text.x=element_text(size=10,colour = "black",angle=90, vjust=0.5),
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
      
      o <- ggplot(data=muts_basis_melt_summary, aes(x=targetfeature, y=percentage,fill=targetstrand))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab(feature)+ylab("Percentage")
      o <- o+geom_errorbar(aes(ymin=percentage-percentage_sd,ymax=percentage+percentage_sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)+scale_y_continuous(labels = scales::percent) #limits=c(0,0.3),breaks=(seq(0,1,0.3)),
      o <- o+scale_x_discrete(limits = as.character(muts_basis_melt_summary$targetfeature))+ggtitle(paste0(outputname,"_signature"))
      #  o <- o+scale_fill_manual(values=mypalette)
      o <- o+theme(axis.text.x=element_text(size=10,colour = "black",angle=90, vjust=0.5),
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
      
      grid.arrange(p,o,ncol=1)
      dev.off()
      
    }
       
  }
}
# Calculate the pvalue from t-test for strand bias (not successful)
denovo_mutscount_Ttest_topography_strandbias <- function(denovo_muts_child_feature, feature,featurelength){
  featurelength$all <-rowSums(featurelength[,2:5])
  
  denovo_muts_child_feature <- denovo_muts_child_feature[denovo_muts_child_feature[,feature] != "others",]
  denovo_muts_child_dis <- data.frame(table(denovo_muts_child_feature[,"Sample"],denovo_muts_child_feature[,feature]))
  names(denovo_muts_child_dis) <- c("Sample",feature, "Freq")
  denovo_muts_child_dis$Sample.Name <- sub("\\_.*","",denovo_muts_child_dis$Sample)
  denovo_muts_child_dis <- merge(denovo_muts_child_dis,samples_details,by="Sample.Name")
  denovo_muts_child_dis$Treatment <- paste0(denovo_muts_child_dis$Sample.Name,"_",denovo_muts_child_dis$Compound.Abbreviation,"_",denovo_muts_child_dis$Concentration)
  denovo_muts_child_dis$targetfeature <- denovo_muts_child_dis[,feature]
  # Control
  control_muts <- denovo_muts_child_dis[denovo_muts_child_dis$Group=="Control",]
  control_profile <- dcast(control_muts[,c("Sample","targetfeature","Freq")], targetfeature~Sample)
  rownames(control_profile) <- control_profile$targetfeature
  
  
  # Compounds
  allcompound_profile <- denovo_muts_child_dis[denovo_muts_child_dis$Group!="Control",]
  compoundlist <- data.frame(table(allcompound_profile$Sample.Name))
  names(compoundlist) <- c("Sample.Name","Freq")
  #compoundlist$P_occurance <- 0
  mut_template <- data.frame(control_profile$targetfeature)
  names(mut_template) <- "targetfeature"
  mutation_type_all <- NULL
  for(i in 1:dim(compoundlist)[1]){
    print(i)
    currentcompound <- denovo_muts_child_dis[denovo_muts_child_dis$Sample.Name==compoundlist[i,"Sample.Name"],]
    compound_profile <- dcast(currentcompound[,c("Sample","targetfeature","Freq")], targetfeature~Sample)
    rownames(compound_profile) <- compound_profile$targetfeature
    
    # calculate pvalue
    compound_bootstrapexposure <- as.data.frame(RemoveBackground_exposure(control_profile[,-1],compound_profile[,-1],1000,1.65))
    colnames(compound_bootstrapexposure) <- paste0("e",seq(1,dim(compound_bootstrapexposure)[2]))
    compound_bootstrapexposure$Mutation <- substr(rownames(compound_bootstrapexposure),1,3)
    compound_bootstrapexposure$Strand <- substr(rownames(compound_bootstrapexposure),5,6)
    compound_bootstrapexposure_melt <- melt(compound_bootstrapexposure,id=c("Mutation","Strand"))
    names(compound_bootstrapexposure_melt) <- c("Mutation","Strand","Sample","exposure")
    chosen_compound_dcast <- dcast(compound_bootstrapexposure_melt,Sample+Mutation~Strand,value.var="exposure")
    
    
    mutation_type <- data.frame(table(chosen_compound_dcast$Mutation))
    names(mutation_type) <- c("Mutation","freq")
    mutation_type$Sample.Name <- as.character(compoundlist[i,"Sample.Name"])
    mutation_type$ratio_mean <- 0.5
    mutation_type$t_pvalue <- 1
    for(j in 1:dim(mutation_type)[1]){
      transtrand <- round(chosen_compound_dcast[chosen_compound_dcast$Mutation==as.character(mutation_type[j,"Mutation"]),"-1"],6)
      nontranstrand <- round(chosen_compound_dcast[chosen_compound_dcast$Mutation==as.character(mutation_type[j,"Mutation"]),"1"],6)
      print(max(mean(transtrand),mean(nontranstrand)))
      mutation_type[j,"max_strand"] <- max(mean(transtrand),mean(nontranstrand))
      mutation_type[j,"t_pvalue"] <- t.test(transtrand,nontranstrand,paired=TRUE)$p.value
      mutation_type[j,"ratio_mean"]  <- sum(transtrand)/(sum(nontranstrand)+sum(transtrand))
    }
    mutation_type_all <- rbind(mutation_type_all,mutation_type)
    
  }
  mutation_type_all$t_qvalue <- p.adjust(mutation_type_all$t_pvalue,method = "BH")
  mutation_type_all$t_pvalue <- round(mutation_type_all$t_pvalue,6)
  mutation_type_all$t_qvalue <- round(mutation_type_all$t_qvalue,6)
  
  mutation_type_all$flag <- ""
  mutation_type_all[which(mutation_type_all$t_qvalue<=0.001),"flag"] <- "***"
  mutation_type_all[which(mutation_type_all$t_qvalue<=0.01 & mutation_type_all$t_qvalue > 0.001),"flag"] <- "**"
  mutation_type_all[which(mutation_type_all$t_qvalue<=0.05 & mutation_type_all$t_qvalue > 0.01),"flag"] <- "*"
  
  write.table(mutation_type_all,"Transtrand_6muttype_53treatment_control_bootstrap_ttest.txt",sep = "\t",col.names = T, row.names = F, quote = F)
  
  # plot
  mutation_type_all_info <- merge(mutation_type_all,samples_details,by="Sample.Name")
  mutation_type_all_info$treatment <- paste0(mutation_type_all_info$Group,"_",mutation_type_all_info$Compound.Abbreviation,"_",mutation_type_all_info$Concentration,"_",mutation_type_all_info$Sample.Name)
  mutation_type_all_info <- mutation_type_all_info[order(mutation_type_all_info$Group,mutation_type_all_info$Compound.Abbreviation),]
  filename=paste0("Transtrand_6muttype_53treatment_control_bootstrap_ttest",".pdf")
  pdf(file=filename, onefile=TRUE,width = 8,height = 10)
  g1 <-ggplot(mutation_type_all_info, aes(x=Mutation, y=treatment)) + geom_tile(data = subset(mutation_type_all_info,t_qvalue<0.1&max_strand>=10),aes(fill=ratio_mean),colour="black",size=0.25)
  g1 <- g1+ geom_tile(data = subset(mutation_type_all_info,t_qvalue>=0.1 | max_strand<10),aes(fill=ratio_mean),alpha=0,colour="black",size=0.25)+geom_text(aes(label=paste(flag)),size=5)
  #  g1 <- g1+scale_y_discrete(limits=as.character(mutation_type_all_info$treatment))
  g1 <-g1 +scale_fill_gradient2(limits=c(0.2, 0.8),midpoint = 0.5)
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
  
}  

# Use ExactSig function, count for replication-timing region
denovo_mutscount_sig_topography2 <- function(denovo_muts_child_feature, feature,featurelength,h,w){
  featurelength$all <-rowSums(featurelength[,2:5])
  
  denovo_muts_child_feature <- denovo_muts_child_feature[denovo_muts_child_feature[,feature] != "others",]
  denovo_muts_child_dis <- data.frame(table(denovo_muts_child_feature[,"Sample"],denovo_muts_child_feature[,feature]))
  names(denovo_muts_child_dis) <- c("Sample",feature, "Freq")
  denovo_muts_child_dis$Sample.Name <- sub("\\_.*","",denovo_muts_child_dis$Sample)
  denovo_muts_child_dis <- merge(denovo_muts_child_dis,samples_details,by="Sample.Name")
  #denovo_muts_child_dis$Treatment <- paste0(denovo_muts_child_dis$Sample.Name,"_",denovo_muts_child_dis$Compound.Abbreviation,"_",denovo_muts_child_dis$Concentration)
  denovo_muts_child_dis$targetfeature <- denovo_muts_child_dis[,feature]
  # Control
  control_muts <- denovo_muts_child_dis[denovo_muts_child_dis$Group=="a_Control",]
  control_profile <- dcast(control_muts[,c("Sample","targetfeature","Freq")], targetfeature~Sample,value.var="Freq")
  rownames(control_profile) <- control_profile$targetfeature
  
  
  # Compounds
  allcompound_profile <- denovo_muts_child_dis[denovo_muts_child_dis$Group!="a_Control",]
  compoundlist <- data.frame(table(allcompound_profile$Sample.Name))
  names(compoundlist) <- c("Sample.Name","Freq")
  #compoundlist$P_occurance <- 0
  mut_template <- data.frame(control_profile$targetfeature)
  names(mut_template) <- "targetfeature"
  for(i in 1:dim(compoundlist)[1]){
    print(i)
    currentcompound <- denovo_muts_child_dis[denovo_muts_child_dis$Sample.Name==compoundlist[i,"Sample.Name"],]
    compound_profile <- dcast(currentcompound[,c("Sample","targetfeature","Freq")], targetfeature~Sample,value.var="Freq")
    rownames(compound_profile) <- compound_profile$targetfeature
    #  ExtractSig(control_profile[,-1],compound_profile[,-1],1000,1.65,mut_template,feature,4.2,7,paste0(feature,as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Compound.Abbreviation"]))
    compound_sig <- RemoveBackground_centroid(control_profile[,-1],compound_profile[,-1],1000,1.65)
    compound_sig$targetfeature <- rownames(compound_sig)
    muts_basis_melt_summary <- merge(mut_template, compound_sig,by="targetfeature",all.x=T)
    muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
    muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$targetfeature),]
    
    muts_basis_melt_summary$percentage <- muts_basis_melt_summary[,"centroid"]/sum(muts_basis_melt_summary[,"centroid"])
    muts_basis_melt_summary$percentage_sd <- muts_basis_melt_summary[,"sd"]/sum(muts_basis_melt_summary[,"centroid"])
    
    outputname=paste0(feature,"_",as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Compound.Abbreviation"])
    
    write.table(muts_basis_melt_summary,paste0(outputname, "_exposure.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
    
    #mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#808080")
    
    filename <- paste0(outputname, "_Signature.pdf")
    pdf(file=filename, onefile=TRUE,width=w,height=h)
    p <- ggplot(data=muts_basis_melt_summary, aes(x=targetfeature, y=centroid))+ geom_bar(stat="identity",position="dodge", width=.8,fill="#3cb44b")+xlab(feature)+ylab("Count")
    p <- p+geom_errorbar(aes(ymin=centroid-sd,ymax=centroid+sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)#+scale_y_continuous(limits=c(-1,1),breaks=(seq(-1,1,0.2)),labels = scales::percent)
    p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$targetfeature))+ggtitle(paste0(outputname,"_exposure"))
    #   p <- p+scale_fill_manual(values=mypalette)
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
    
    
    o <- ggplot(data=muts_basis_melt_summary, aes(x=targetfeature, y=percentage))+ geom_bar(stat="identity",position="dodge", width=.8,fill="#3cb44b")+xlab(feature)+ylab("Percentage")
    o <- o+geom_errorbar(aes(ymin=percentage-percentage_sd,ymax=percentage+percentage_sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)+scale_y_continuous(labels = scales::percent) #limits=c(0,0.3),breaks=(seq(0,1,0.3)),
    o <- o+scale_x_discrete(limits = as.character(muts_basis_melt_summary$targetfeature))+ggtitle(paste0(outputname,"_signature"))
    #  o <- o+scale_fill_manual(values=mypalette)
    o <- o+theme(axis.text.x=element_text(size=10,colour = "black"),
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
    
    grid.arrange(p,o,ncol=1)
    dev.off()
    
  }
}  
denovo_6muttype_topography_byTreatment_control_aggragate <- function(denovo_muts_child_feature, feature,featurelength,Isprint,h,w, outputname){
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  featurelength$all <-rowSums(featurelength[,2:5])
  
  denovo_muts_child_feature <- denovo_muts_child_feature[denovo_muts_child_feature[,feature] != "others",]
  denovo_muts_child_feature$Mutation <- paste0(denovo_muts_child_feature$Ref2,">",denovo_muts_child_feature$Alt2)
  denovo_muts_child_dis <- data.frame(table(denovo_muts_child_feature[,"Sample"],denovo_muts_child_feature[,feature], denovo_muts_child_feature[,"Mutation"]))
  names(denovo_muts_child_dis) <- c("Sample",feature,"Mutation", "Freq")
  
  denovo_muts_child_dis$Sample.Name <- "Control"
  
  gtc <- ddply(denovo_muts_child_dis, c(feature,"Mutation","Sample.Name"),summarise, N=length(Sample.Name),sum=sum(Freq))
  gtc$targetfeature <- gtc[,feature]
  featurelist <- data.frame(table(gtc[,feature]))[,1]
  if(Isprint=="Print"){
    filename=paste0(outputname,".pdf")
    pdf(file=filename, onefile=TRUE,height=h,width = w) 
    d1 <- ggplot(gtc,aes(x=Mutation,y=sum,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
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
    d1 <- d1+facet_wrap(~Sample.Name,ncol=5)
    print(d1)
    dev.off()
    
    
  }
  write.table(gtc,paste0(outputname,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  return(gtc)
}



