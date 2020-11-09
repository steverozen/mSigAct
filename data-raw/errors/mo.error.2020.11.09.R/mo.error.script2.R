setwd("~/mSigAct/data-raw/errors/mo.error.2020.11.09.R")
library(mSigAct)

sigminer.optimal.exp.DBS.example <- ICAMSxtra::ReadExposure("sigminer.assignment.csv")
hdp.DBS.signatures <- ICAMS::ReadCatalog("hdp.DBS.signatures.csv")
monster.DBS.catalog.example <- ICAMS::ReadCatalog("monster.DBS.example.catalog.csv")


mSigAct.optimal.exp.DBS.example <- sigminer.optimal.exp.DBS.example
mSigAct.optimal.exp.DBS.example[,1:ncol(mSigAct.optimal.exp.DBS.example)] <- 0

m.opts <- DefaultManyOpts()
sample <- "Adrenal::POG570_P18843" #the sample can cause error

#select signatures present in sigminer assignment
selected.sigs <- row.names(sigminer.optimal.exp.DBS.example)[(which(sigminer.optimal.exp.DBS.example[,sample]>0))]

sp <- monster.DBS.catalog.example[,sample,drop=F]
sg <- hdp.DBS.signatures[,selected.sigs,drop=F]

SA.out <-
  SparseAssignActivity(spectra              = sp,
                       sigs                 = sg,
                       # eval_f               = ObjFnBinomMaxLHNoRoundOK,
                       eval_f               = ObjFnBinomMaxLHMustRound,
                       m.opts               = m.opts,
                       mc.cores.per.sample  = 1,
                       num.parallel.samples = 1,
                       max.level            = 5, # I think this is the default
                       p.thresh              = 0.05 # You can test different values to get more ore less sparsity; lower values will give more sparsity
  )
mSigAct.optimal.exp.DBS.example[,sample] <- 0
mSigAct.optimal.exp.DBS.example[selected.sigs,sample] <- SA.out$exposure


