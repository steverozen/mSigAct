## Supplementary code 1, Boot et al., Script 2

## Mutational signature assignment on SBS mutation spectra of 
## Duplex Sequencing data of UTUCs, normal tissues, and DNA from blood and urine

recompute.all <- FALSE

########################################################################## 
################ load packages, helper functions and data ################
########################################################################## 

source("Helper functions.R")
require.all.packages()

# load the mutation spectra of our samples
cat<-ICAMS::ReadCatalog("input/Twinstrand_Catalog_96_noGermline.csv")


########################################################################## 
####################### load reference signatures ########################
########################################################################## 

# load reference signatures (COSMIC v3.1)
# Source: https://cancer.sanger.ac.uk/cosmic/signatures/SBS/index.tt
ref.sigs <- ICAMS::ReadCatalog("input/COSMIC_Mutational_Signatures_v3.1_SBS.csv")
ref.sigs <- ICAMS::as.catalog(ref.sigs,
                              ref.genome = "GRCh37",
                              region = "genome",
                              catalog.type = "counts.signature")

# The Duplex Sequencing only targets a small part of the human genome
# calculate the trinucleotide composition of the DS target region
# and adjust the reference signatures to take this difference into account
DS.target.BED<-ICAMS:::ReadBedRanges("input/Twinstrand.bed")
DS_abundance<-ICAMS:::GetCustomKmerCounts(k = 3,
                                          ref.genome = "GRCh38",
                                          custom.ranges = DS.target.BED)
DS_abundance<-ICAMS:::Collapse192AbundanceTo96(DS_abundance[,1])
ref.sigs.DS<-ICAMS::TransformCatalog(ref.sigs,
                                     target.ref.genome = NULL,
                                     target.region = "unknown",
                                     target.catalog.type = "counts.signature",
                                     target.abundance = DS_abundance)

ICAMS::WriteCatalog(ref.sigs.DS,file="input/ref.sigs.DS.csv")

# We use all SBS signatures assigned to the WGS data as prior (read expected) signatures
# The "output/SBS_assignment.RData" was generated in 1) WGS_analysis.R
load("output_MAP/SBS_assignment_WGS.RData")
prior.sigs <- rownames(MAP.out$proposed.assignment)[rowSums(MAP.out$proposed.assignment) > 0]

########################################################################## 
################################ Run HDP ################################# 
########################################################################## 

## only run if required
if (recompute.all || !file.exists("HDP_DS/hdp.retval.Rdata")) {
  message("Running RunHdpParallel")
  
  seed <- 1324  
  num.jobs <- 20
  
  mSigHdp::RunHdpxParallel(
    input.catalog       = "input/Twinstrand_Catalog_96_noGermline.csv",
    seedNumber          = seed,
    K.guess             = 20 ,
    multi.types         = TRUE,
    verbose             = TRUE,
### post.burnin         = 15000,
    post.n              = 200, 
    post.space          = 300,
    post.cpiter         = 3,
    post.verbosity      = 0 ,
    CPU.cores           = 20, 
    num.child.process   = 20,
    cos.merge           = 0.9,
### min.sample          = 1,
    ground.truth.sig    = "input/ref.sigs.DS.csv",
    ground.truth.exp    = NULL,
    overwrite           = TRUE,
    out.dir             = "HDP_DS_MAP",
    gamma.alpha         = 1,
    gamma.beta          = 50,
    prior.sigs          = ref.sigs[,prior.sigs],
    prior.pseudoc       = round(rep(median(colSums(cat)),length(prior.sigs)),0))
}

# We extracted one signature that is ubiquitous to the DS data: hdp.N1
# This represents either a true mutational process or a DS specific artifact.
# Either way, we will use it as an input signature for the signature assignment


########################################################################## 
###################### Define Signature universe ######################### 
########################################################################## 

# add the DS-specific mutational signature for signature assignment
# hdp.2 is the DS-specific signature
HDPsigs<-ICAMS::ReadCatalog(paste0("HDP_DS_MAP/extracted.signatures.csv"))
HDPsigs<-ICAMS::as.catalog(HDPsigs,
                           ref.genome = "GRCh37",
                           region = "genome",
                           catalog.type = "counts.signature")

# Using all DBS signatures for reconstruction of DBS mutation spectra
# This allows all previously described DBS signatures to compete with
# the newly discovered signature (hdp.2)
ref.sigs.DS<-cbind(ref.sigs.DS,HDPsigs[,"hdp.2",drop=F])

# The blood sample of CGMH-AB18 shows obvious cisplatin mutagenesis.
# Therefore we add SBS31 in the analysis
UTUC_universe<-c(prior.sigs,"hdp.2","SBS31")


########################################################################## 
################# perform mutational signature assignment ################
########################################################################## 

if(recompute.all || !file.exists("output_MAP/SBS_assignment_DS.RData")){
  
  ## make up sigs.prop
  sigs.prop<-rep(0.5,length(UTUC_universe))
  names(sigs.prop)<-UTUC_universe
    
  MAP.out<- 
    MAPAssignActivity(spectra = cat,
                      sigs = ref.sigs.DS[,UTUC_universe],
                      sigs.presence.prop = sigs.prop,
                      p.thresh = 0.05/length(sigs.prop),
                      max.level = length(sigs.prop)-1,
                      num.parallel.samples = 1,
                      mc.cores.per.sample = 20,
                      output.dir = paste0(getwd(),"/output_MAP/"))
  
  # Save the signature assignments
  save(MAP.out,file="output_MAP/SBS_assignment_DS.RData")
} else {
  load("output_MAP/SBS_assignment_DS.RData")
}



########################################################################## 
########################## plot the assignments ########################## 
########################################################################## 

sparseAssign<-sparseAssign[c("hdp.2",paste0("SBS",c(1,5,40,2,13,8,31,22))),]

cairo_pdf("output_MAP/SBS_assignment_DS.pdf", width=7, height=8, onefile=T)
par(mfrow=c(3,1), oma=c(2,2,2,2), mar=c(5,4,4,1),lwd = 0.3)

# start with plotting the total mutation load
mutLoads<-colSums(sparseAssign)
names(mutLoads)<-colnames(sparseAssign)

# plot aesthetics
cols<-rep("black",length(mutLoads))
cols[substr(colnames(sparseAssign),1,1) == "N"]<-"grey50"
cols[substr(colnames(sparseAssign),1,1) == "U"]<-"orange"
cols[substr(colnames(sparseAssign),1,1) == "B"]<-"red"
xlimXtra = 1.5
xlim=xlim=c(0.5, length(mutLoads)*xlimXtra)

bp<-barplot(mutLoads,log="y",xlim=xlim,col=cols,
            xlab="",ylab="Total SBS load",
            main="SBS mutation load WGS data",
            las=1,xaxt="n",yaxt="n",ylim=c(10,2500))
addXaxis(bp,labels=reconstructionAccuracy[,1])
addYaxis(c(10,25,50,100,250,500,1000,2500))
legend(ncol(sparseAssign)*1.2-0.25,max(mutLoads,na.rm = T),
       legend=c("Tumor","Normal","Urine","Blood"),
       ncol=1, fill=c("black","grey50","orange","red"),
       cex = 1, bty='n')

sparseAssign2<-sweep(sparseAssign,2,colSums(sparseAssign),"/")

plot.stacked(as.matrix(sparseAssign2),
             main = "Signature assignments WGS samples",
             cols=c("lemonchiffon","grey40","grey60","grey80",
                    "#009900","#4bd24b",
                    "orange","blue",
                    "red"),
             ylab="Proportion",
             las=2,cex.axis=1,cex.names=1,
             xlimXtra = xlimXtra)

# Plot the reconstruction accuracy
plot(bp,reconstructionAccuracy[,2],col=cols,pch=16,
     xlim=xlim,ylim=c(0,1),xaxt="n",yaxt="n",bty='n',
     xlab="",ylab="reconstruction accuracy (cosine)",
     main="Reconstruction accuracy of SBS spectra")
addXaxis(bp,labels=reconstructionAccuracy[,1])
addYaxis(seq(0,1,0.2))
legend(ncol(sparseAssign)*1.2-0.25,1,ncol=1, 
       legend=c("Tumor","Normal","Urine","Blood"),
       col=c("black","grey50","orange","red"), pch=c(16,16),
       cex = 1, bty='n')

dev.off()

############################# End of script ############################## 
########################################################################## 