MSISigs <- function() {
  SBS <- paste0("SBS", c(6, 14, 15, 20, 21, 26, 44))
  DBS <- paste0("DBS", c(7, 10))
  ID <- paste0("ID", c(1, 2, 7))
  return(list(SBS = SBS, DBS = DBS, ID = ID))
}