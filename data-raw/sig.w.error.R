## stripped down script to plot signatures with error-bars

## load the helper functions
source("Helper_functions_signature_stability.R")

## NB: the helper function is the ICAMS signature plot function with 2 modifications:
##          - it returns the output from the bp call used to plot the signature
##          - it has an additional input parameter 'extraY' which allows me to adjust the Y-axis to make room for the error bar

## load catalog
load("background_spectra.Rdata")

## convert to proportions and calculate mean + sd
catPerc<-as.data.frame(sweep(catSBS$catSBS96,2,colSums(catSBS$catSBS96),"/"))
catPerc$HepG2_avg<-apply(catPerc[,1:3],1,FUN = mean)
catPerc$HepG2_sd<-apply(catPerc[,1:3],1,FUN = sd)

## convert to ICAMS catalog object
catPerc2<-ICAMS::as.catalog(catPerc,"GRCh37","genome","counts.signature")

## plot the signature (the average of the clones)
tmp<-catPerc2[,7,drop=F]
i<-my.PlotCatalog.SBS96Catalog(tmp,extraY=1.05)

## NB: i contains the output from the bp function used to plot the signature
## this output is used in the next step for the x-coordinates of the error bars
## inside the line below, the sd is converted to SEM, by dividing by sqrt(3)

suppressWarnings(arrows(x0=i$plot.object,
                        y0=catPerc2[,7]+(catPerc2[,8]/sqrt(3)), # Up
                        y1=catPerc2[,7]-(catPerc2[,8]/sqrt(3)), # down
                        angle=90, # flat
                        code=3,
                        length=0.025 #length of tick
                        ))

