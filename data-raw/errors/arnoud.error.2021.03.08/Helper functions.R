## Helper functions for Supplementary code 1 to Boot et al.

## Require CRAN library and if not installed prompt for installation
check.install.CRAN <- function(pp) {
  if (!suppressWarnings(require(pp, character.only = TRUE, quietly = TRUE))) {
    message("package ", pp, " is not installed")
    message("to install from CRAN do:\ninstall.packages(\"", pp, "\")\n")
    return(FALSE)
  }
  return(TRUE)
}

## Require github library and if not installed prompt for installation
check.install.github <- function(pp.dot.ref) {
  pp <- pp.dot.ref$repo
  if (!suppressWarnings(require(pp, character.only = TRUE, quietly = TRUE))) {
    message("\npackage ", pp, " is not installed")
    message("to install from github do:")
    message("remotes::install_github(\"steverozen/", pp,
            "\", ref = \"", pp.dot.ref$ref, "\")\n")
    return(FALSE)
  }
  return(TRUE)
}

require.all.packages <- function() {
  r1 <- unlist(lapply(c("lsa", "RColorBrewer", "remotes"), check.install.CRAN))
  r2 <- unlist(
    lapply(
      list(
        # Note: SynSigEval is suggested by mSigHdp, but its functionality is not
        # needed for these scripts
        list(repo = "mSigHdp", ref = "master"),
        list(repo = "mSigBG",  ref = "master"),
        list(repo = "mSigAct", ref = "master"),
        list(repo = "ICAMS",   ref = "master")),
      check.install.github
    ))
  if (!all(c(r1, r2))) stop("Please install necessary packages")
}

## plotting mutational signature assignments
plot.stacked <- function(proportions, main,leg.adj=1,cols=NULL,
                         xlimXtra=1,density=NULL,angle=NULL,...) {
  if(is.null(cols)){
    require(RColorBrewer)
    palette <- brewer.pal(12, 'Paired')
    palette <- rep(palette,3)
    palette <- rev(palette[1:nrow(proportions)])
    palette[palette == "#E31A1C"]<-"grey50"
    palette[length(palette)]<-"red"
  } else {
    palette<-cols
  }
  
  bp <- barplot(
    proportions,
    main=main,
    xlim = c(1, ncol(proportions)*xlimXtra), 
    col = palette,
    border = "black",
    las = 3,
    cex.main=0.8,
    angle=angle,
    density=density,
    ...)
  legend(ncol(proportions)*1.2-0.25*leg.adj,1,
         legend=rev(rownames(proportions)),
         angle=rev(angle), 
         density=rev(density),
         ncol=1, fill=rev(palette),
         cex = 1, bty='n')
}

recAccuracy<-function(assignments=NULL,cat=NULL,sigs=NULL){

  if(is.null(assignments)|is.null(cat)|is.null(sigs)){
    stop("Assignments, cat and sigs must all be provided")
  }
  
  if(!sum(rownames(assignments) %in% colnames(sigs)) == nrow(assignments)){
    stop("Not all signatures in assignments are present in the sigs object")
  }
  
  if(!sum(colnames(cat) %in% colnames(assignments)) == ncol(cat)){
    stop("colnames cat and assignments must be identical")
  }
  
  output<-NULL
  for(smp in colnames(cat)){
    x<-cat[,smp]
    
    y<-rep(0,nrow(cat))
    for(i in which(assignments[,smp] > 0)){
      y<-y+assignments[i,smp]*sigs[,rownames(assignments)[i]]
    }
    
    z<-round(as.numeric(cosine(x,y)),3)
    
    output<-rbind(output,c(sample=smp,reconstructionAccuracy=z))
  }
  return(output)
}

addXaxis<-function(bp,labels,...){
  axis(1,at=bp,tick = T,
       labels =labels,las=2,cex=0.8,...)
  axis(1,at=c(-1,bp,max(bp)+1),tick = T,labels =F,tck=-0.001,...)
}

addYaxis<-function(labels,s=1,...){
  axis(2,at=labels*s,tick = T,
       labels =labels,las=1,cex=0.8,...)
}

