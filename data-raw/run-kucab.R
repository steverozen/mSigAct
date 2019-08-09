

setwd("data-raw/Cell_MutagenSig-master/")
source("Figure3_SubstitutionSig.R")

# Variable control_profile contains the control spectra

Unstaple96 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           substr(c1, 3, 3),
           substr(c1, 7, 7),
           substr(c1, 5, 5))
  return(retval)
}

Restaple96 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           "[",
           substr(c1, 2, 2),
           ">",
           substr(c1, 4, 4),
           "]",
           substr(c1, 3, 3))
  return(retval)
}
           
kucab.controls <- control_profile  
rownames(kucab.controls) <- Unstaple96(rownames(kucab.controls))
kucab.controls <- 
  ICAMS::as.catalog(kucab.controls, catalog.type = "counts", region = "genome")
control.details <- 
  muts_summary_details[muts_summary_details$Group == "a_Control", ]
new.colnames <- 
    cbind(colnames(kucab.controls),
    control.details[colnames(kucab.controls),
                    c("Compound", "Concentration")])
colnames(kucab.controls) <- 
  apply(new.colnames, MARGIN = 1, FUN = paste0, collapse = "_")
kucab.controls <- kucab.controls[ICAMS::catalog.row.order$SBS96, ]
usethis::use_data(kucab.controls, overwrite = TRUE)
ICAMS::PlotCatalogToPdf(mSigAct::kucab.controls, "kucab.controls.pdf")

kucab.compute.control.dist <- function(num.controls) {
  x <- colSums(mSigAct::kucab.controls)
  retval <- NULL
  set.seed(217)
  for(j in 1:10000){
    selx <- sample(x, num.controls, replace = TRUE)
    retval <- c(retval, mean(selx))
  }
  return(retval)
}

kucab.control.dist <- list()
for (i in 2:4) kucab.control.dist[[i]] <- kucab.compute.control.dist(i)
usethis::use_data(kucab.control.dist)
