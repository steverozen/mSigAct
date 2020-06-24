dir.create("Figure5_results")
setwd("./Figure5_results/")
source("../Header.R")
samples_details <- read.table("../00_data/final_mutagen_info_forR_v4_u.txt",sep = "\t",header = T,as.is = T, quote="\"")

indel.classified <- read.table("../00_data/denovo_subclone_indels.final.txt",sep = "\t",header = T, as.is = T)
indel.classified$Sample.Name <- sub("\\_.*","",indel.classified$Sample)
indel.classified <- indel.classified[indel.classified$Sample.Name!="MSM0",]
indel.classified_details <- merge(indel.classified, samples_details, by="Sample.Name")
indel.classified_details <- indel.classified_details[indel.classified_details$VAF.Tum_Cal>=0.2,]
muts_template <- read.table("../00_data/indel_29channels_template.txt",sep = "\t",header = T, as.is = T)

###############################################
#
# Indel catalogue 29 channels 

#[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
#[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
#[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
#[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
#[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
#[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh

###############################################
indel.classified_details2 <- indel.classified_details
#indel.classified_details2$length_type <- indel.classified_details2$indel.length
#indel.classified_details2[indel.classified_details2$indel.length>1,"length_type"] <- ">1bp"
#indel.classified_details2$Subtype <- paste0(indel.classified_details2$Type,"_",indel.classified_details2$classification,"_",indel.classified_details2$length_type)
#indel.classified_details2[indel.classified_details2$indel.length==1,"Subtype"] <- paste0(indel.classified_details2[indel.classified_details2$indel.length==1,"Subtype"], "_",indel.classified_details2[indel.classified_details2$indel.length==1,"change.pyr"])
indel.classified_details2$Subtype <- NULL
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[+C]A"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[+C]G"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="T","Subtype"] <- "[+C]T"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==1,"Subtype"] <- "[+C]C"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==2,"Subtype"] <- "[+C]CC"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount>2,"Subtype"] <- "[+C]LongRep"

indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[+T]A"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="C","Subtype"] <- "[+T]C"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[+T]G"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==1,"Subtype"] <- "[+T]T"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==2,"Subtype"] <- "[+T]TT"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount>2,"Subtype"] <- "[+T]LongRep"

indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount<1,"Subtype"] <- "[+>1]NonRep"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount>=1,"Subtype"] <- "[+>1]Rep"

indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[-C]A"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[-C]G"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="T","Subtype"] <- "[-C]T"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==1,"Subtype"] <- "[-C]C"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==2,"Subtype"] <- "[-C]CC"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount>2,"Subtype"] <- "[-C]LongRep"


#indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification=="Repeat-mediated" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount>1,"Subtype"] <- "[-C]_Rep"

indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[-T]A"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="C","Subtype"] <- "[-T]C"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[-T]G"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==1,"Subtype"] <- "[-T]T"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==2,"Subtype"] <- "[-T]TT"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount>2,"Subtype"] <- "[-T]LongRep"
#indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification=="Repeat-mediated" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount>1,"Subtype"] <- "[-T]_Rep"

indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification !="Microhomology-mediated" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount<1,"Subtype"] <- "[->1]Others"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification=="Repeat-mediated" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount>=1,"Subtype"] <- "[->1]Rep"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification =="Microhomology-mediated","Subtype"] <- "[->1]Mh"



indel_catalogue <- data.frame(table(indel.classified_details2$Sample,indel.classified_details2$Subtype))
names(indel_catalogue) <- c("suclone","subtype","freq")
indel_catalogue <- dcast(indel_catalogue,suclone~subtype)
indel_catalogue$Sample.Name <- sub("\\_.*","",indel_catalogue$suclone)
write.table(indel_catalogue, paste0("indel_catalogue",".txt"),sep = "\t",col.names = T, row.names = F, quote = T)

indel_catalogue <- merge(indel_catalogue,samples_details, by="Sample.Name", all.x=T)
indel_catalogue$indel_num <- rowSums(indel_catalogue[,3:(2+29)])
write.table(indel_catalogue, paste0("indel_catalogue_withinfo",".txt"),sep = "\t",col.names = T, row.names = F, quote = F)

muts_control <- indel_catalogue[indel_catalogue$Group=="a_Control",]
muts_compound <- indel_catalogue[!indel_catalogue$Group=="a_Control",]



##########################
# Figure 5A
##########################
indel_catalogue <- read.table("indel_catalogue_withinfo.txt",sep = "\t",header = T, as.is = T, quote = "\"",check.names = F)
muts_control <- indel_catalogue[indel_catalogue$Group=="a_Control",]
muts_compound <- indel_catalogue[indel_catalogue$Group!="a_Control",]

# Control
parentmuts <- t(muts_control[,3:(2+29)])
colnames(parentmuts) <- muts_control$suclone
parentmuts <- as.data.frame(parentmuts)

parentmuts$indelsubtype <- rownames(parentmuts)
control_mean <- plotCountbasis_aggragated_indel_29types_6(parentmuts,3,10,"Fig5A")
write.table(control_mean, paste0("Fig5A",".txt"),sep = "\t",col.names = T, row.names = F, quote = T)

##########################
# Figure 5B
##########################
chosen_control_all <- NULL
for(j in 1:1000){
  chosen_control <- apply(muts_control[,3:(2+29)],2,function(x) x[sample(1:dim(muts_control[,3:(2+29)])[1],1,replace=T)])
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
  currentcompound <- muts_compound[muts_compound$Sample.Name==as.character(compoundlist[i,1]),]
  
  # permutation
  chosen_compound_all <- NULL
  for(j in 1:1000){
    chosen_compound <- apply(currentcompound[,3:(2+29)],2,function(x) x[sample(1:dim(currentcompound[,3:(2+29)])[1],1,replace=T)])
    chosen_compound_all <- rbind(chosen_compound_all,chosen_compound)
  }
  
  currentcompound_muts <- t(chosen_compound_all)
  currentcompound_muts <- as.data.frame(currentcompound_muts)
  
  centroid_compound <- rowMeans(currentcompound_muts)
  sd_compound <- sd_highD(currentcompound_muts,centroid_compound-centroid_control)
  sd_control<- sd_highD(chosen_control_muts,centroid_compound-centroid_control)
  
  compoundlist[i,"profile_SNR"]  <- norm(as.matrix(centroid_compound-centroid_control),"f")/(sd_compound+sd_control)
}
Sample_withSig_SNR <- read.table("../Figure2_results/Fig2C.txt",sep = "\t",header = T, as.is = T, quote="\"")
Sample_withSig_SNR2 <- merge(Sample_withSig_SNR,compoundlist[,c("Sample.Name","profile_SNR")],by="Sample.Name")
write.table(Sample_withSig_SNR2,"Sample_withSig_SNR2.txt",sep = "\t",col.names = T, row.names = F, quote = F)



###############################################
# Figure 5B
# Exposure and Signature
###############################################

indel_catalogue <- read.table("./indel_catalogue_withinfo.txt",sep = "\t",header = T, as.is = T, quote = "\"",check.names = F)

muts_control <- indel_catalogue[indel_catalogue$Group=="a_Control",]
muts_compound <- indel_catalogue[!indel_catalogue$Group=="a_Control",]


# control
controlmuts <- t(muts_control[,3:(2+29)])
colnames(controlmuts) <- muts_control$suclone
controlmuts <- as.data.frame(controlmuts)
controlmuts$indelsubtype <- rownames(controlmuts)
controlclones <- controlmuts[,-ncol(controlmuts)]


#############################################
# 1) Extract signature for each compound condition

Sample_withSig <- read.table("./Sample_withSig_SNR2.txt",sep = "\t",header = T, as.is = T,quote = "\"")
sig_indels_full <- Sample_withSig[Sample_withSig$pvalue<=0.01 & Sample_withSig$profile_SNR>=2 & Sample_withSig$mean>=20,]



compoundlist <- data.frame(table(muts_compound$Sample.Name))
names(compoundlist) <- c("Sample.Name","Freq")
compoundlist <- compoundlist[compoundlist$Sample.Name %in% sig_indels_full$Sample.Name,]
stability_all <- NULL
for(i in 1:dim(compoundlist)[1]){
  
  print(i)
  currentcompound <- muts_compound[muts_compound$Sample.Name==as.character(compoundlist[i,1]),]
  
  currentcompound_muts <- currentcompound[,3:(2+29)]
  currentcompound_muts <- t(currentcompound_muts)
  colnames(currentcompound_muts) <- currentcompound$suclone
  currentcompound_muts <- as.data.frame(currentcompound_muts)
  currentcompound_muts$indelsubtype <- rownames(currentcompound_muts)
  compoundclones <- currentcompound_muts[,-ncol(currentcompound_muts)]
  
  currentcompound_stability <- ExtractSig_aggregate_indels(controlclones,compoundclones,1.28,muts_template,5,6,paste0("indel_",as.character(compoundlist[i,1]),"_",samples_details[samples_details$Sample.Name==as.character(compoundlist[i,1]),"Compound.Abbreviation"]))
  stability_all <- rbind(stability_all,c(currentcompound_stability))
  
  
  #MCSignature_indels_2(controlclones,compoundclones,1000,1.65,5,3.3,paste0("indels_",as.character(compoundlist[i,1]),"_",currentcompound[1,"Compound.Abbreviation"]))
  # current_P_occurance <-   MCSignature_v2(controlclones,compoundclones,1000,5,3.3,paste0("indels_",as.character(compoundlist[i,1]),"_",currentcompound[1,"Compound.Abbreviation"]))
  
  # compoundlist[i,"P_occurance"] <- current_P_occurance
  #ChemExposure(controlclones,compoundclones,6,5,paste0("indels_",as.character(compoundlist[i,1]),"_",currentcompound[1,"Compound.Abbreviation"],"_exposure"))
}
#control_profile_set=controlclones
#compound_profile_set=compoundclones
stability_all_2 <- cbind(stability_all,compoundlist)
stability_all_2 <- as.data.frame(stability_all_2)
names(stability_all_2) <- c("min_simi","max_simi","Sample.Name","subclone_num")
Sample_withSig_SNR3 <- merge(Sample_withSig_SNR2,stability_all_2,by="Sample.Name")
Sample_withSig_SNR3$stability <- round(Sample_withSig_SNR3$max_simi,2)
write.table(Sample_withSig_SNR3,"Fig5B.txt",sep = "\t",col.names = T, row.names = F, quote = F)

# Ten of them shown in Figure 5B have SNR >= 2, average indel number per subclone >=20 and stability >=0.7. 

###################################################################
#
# cosine similarities between signatures 
# including lung cancer signature from pancan
#
###################################################################
lung_sig <- read.table("../00_data/lung_signature.tsv",sep = "\t", header = T, as.is = T)
names(lung_sig) <- c("indelsubtype","indeltype","lung_percentage")

plotSinglebasis_indel_29types_6 <- function(muts_basis,h,w,outputname){
  muts_template <- read.table("../00_data/indel_29channels_template.txt",sep = "\t",header = T, as.is = T)
  muts_basis <- merge(muts_template, muts_basis[,c(1,3)],by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
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
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=percentage,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
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
plotSinglebasis_indel_29types_6(lung_sig,2.8,6,"FigS4_lung_sig")

sig_sample <- c("MSM0.110","MSM0.83","MSM0.92","MSM0.2","MSM0.103","MSM0.95","MSM0.87","MSM0.93")
sigs_12_full <- samples_details[samples_details$Sample.Name%in%sig_sample,]
sig_all <- NULL
for(i in 1:length(sig_sample)){
  sig_file <- read.table(paste0("indel_",sigs_12_full[i,"Sample.Name"],"_",sigs_12_full[i,"Compound.Abbreviation"],"_exposure.txt"), sep = "\t", header = T, as.is = T)
  sig_all <- cbind(sig_all,sig_file[,"percentage"])
  
}
sig_all <- as.data.frame(sig_all)

sig_all$indelsubtype <- sig_file[,"indelsubtype"]
sig_all <- merge(sig_all,lung_sig)

colnames(sig_all)[2:9] <- sigs_12_full$Compound.Abbreviation


for(i in 2:9){
  print(cos_similarity(sig_all[,"percentage"],sig_all[,i]))
}
#[1] 0.2452572 DMH
#[1] 0.3649882 Cisplatin
#[1] 0.9177187 BaP
#[1] 0.9279615 BPDE
#[1] 0.8100585 DBADE
#[1] 0.5553006 1,8-DNP
#[1] 0.6245682 3-NBA
#[1] 0.5325716 6-Nitrochrysene


control_sig <- read.table("./Fig5A.txt",sep = "\t", header = T, as.is = T)
control_sig <- control_sig[,c(1,38)]
control_sig$control_percentage <- control_sig[,2]/sum(control_sig[,2])

sig_all <- merge(sig_all,control_sig[,c(1,3)],by="indelsubtype")
sig_all <- sig_all[,c(2:9,11,12)]
sig_all <- as.data.frame(t(sig_all))
cossimil <- as.matrix(proxy::simil(as.matrix(sig_all), diag=FALSE, upper=FALSE, method="cosine",auto_convert_data_frames = FALSE))
cossimil[lower.tri(cossimil,diag=TRUE)]=NA
cossimil <- as.data.frame(as.table(cossimil))
cossimil=na.omit(cossimil)
names(cossimil) <- c("sample1","sample2","simil")
cossimil$simil <- round(cossimil$simil,2)
write.table(cossimil,"Fig5C_cossim_indel.txt",sep="\t",col.names = T, row.names = F, quote = F)
filename=paste0("Fig5C_cossim_indel_red",".pdf")
pdf(file=filename, onefile=TRUE,width = 5,height =4)
g1 <-ggplot(cossimil, aes(x=sample1, y=sample2)) + geom_tile(aes(fill=simil),colour="white")+geom_text(aes(label=paste(simil)),size=3)
g1 <-g1 +scale_fill_gradient2(high="red", low="white",space="Lab",limits=c(0, 1))
#g1 <-g1 +scale_x_discrete(limits = unique(as.character(cossimil$sample1)))
#g1 <-g1 +scale_y_discrete(limits = unique(as.character(cossimil$sample2[order(as.character(cossimil$sample2))])))

g1 <-g1 +theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
g1
dev.off()

##################
# Fig5C Cisplatin
##################
cisplatin_indels <- indel.classified_details[indel.classified_details$Sample.Name=="MSM0.83",]
write.table(cisplatin_indels,"cisplatin_indel_83.txt", sep = "\t", col.names = T, row.names = F, quote = F)
# Using Excel to count and make figure
#           TA	TC	TG	TT	T3	T4	T5	T6	T7	T8	T9
# GG	       1   5	18	9	  10	5	  3	  0	  0	  0 	0
# Others	   0	 1	1	  2	  2	  2	  0	  1	  2	  3	  3

##################
# Fig5D DBADE
##################
indel_mutagen <- indel.classified_details[indel.classified_details$Compound.Abbreviation=="DBADE",]
write.table(indel_mutagen,"DBADE_indel.txt",sep = "\t",col.names = T, row.names = F, quote = F)
# Using Excel to count and make figure
#   AA	AG	GA	GC	GG	G3	G4	G5	Others
#    2	1	   1	1	  11	30	9	  3	  0
#    0	0	   0	0	  0	  8	  0	  1	  0








