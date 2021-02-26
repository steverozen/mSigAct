# Source this file from mSigAct package root

cat(getwd(), "\n")

tmp.SBS96 <- 
  data.table::fread("data-raw/package-data-related/source-files/COSMIC-v3-SBS-proposed-aetiology-short-form.csv")
sigs.etiologies.SBS96 <- 
  matrix(tmp.SBS96$proposed.aetiology, ncol = 1, 
         dimnames = list(tmp.SBS96$name, "proposed.etiology"))

tmp.SBS192 <- 
  data.table::fread("data-raw/package-data-related/source-files/COSMIC-v3-SBS-proposed-aetiology-short-form.csv")
sigs.etiologies.SBS192 <- 
  matrix(tmp.SBS192$proposed.aetiology, ncol = 1, 
         dimnames = list(AddMinusEFor192(tmp.SBS96$name), "proposed.etiology"))

tmp.DBS78 <- 
  data.table::fread("data-raw/package-data-related/source-files/COSMIC-v3-DBS-proposed-aetiology-short-form.csv")
sigs.etiologies.DBS78 <- 
  matrix(tmp.DBS78$proposed.aetiology, ncol = 1, 
         dimnames = list(tmp.DBS78$name, "proposed.etiology"))

tmp.ID <- 
  data.table::fread("data-raw/package-data-related/source-files/COSMIC-v3-ID-proposed-aetiology-short-form.csv")
sigs.etiologies.ID <- 
  matrix(tmp.ID$proposed.aetiology, ncol = 1, 
         dimnames = list(tmp.ID$name, "proposed.etiology"))

sigs.etiologies <- list(SBS96 = sigs.etiologies.SBS96,
                        SBS192 = sigs.etiologies.SBS192,
                        DBS78 = sigs.etiologies.DBS78,
                        ID = sigs.etiologies.ID)

usethis::use_data(sigs.etiologies, internal = TRUE, overwrite = TRUE)
