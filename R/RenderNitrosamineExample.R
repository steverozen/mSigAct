RenderNitroExample <- 
  function(whichnitro, bg.inflate.factor = 1, count.nbinom.size = 20) {
    
    outfile <- file.path(
      devtools::package_file("data-raw"),
      paste0("Nitrosamine-", whichnitro, 
             "-x", bg.inflate.factor, 
             "-disp", count.nbinom.size,
             ".html"))
    
    rmarkdown::render(
      
      input          = devtools::package_file("data-raw/NitrosamineExample.Rmd"),
      
      output_file    = outfile,
      
      output_options = list(title = paste("Background subtraction for", 
                                          whichnitro,
                                          "with bg inflation",
                                          bg.inflate.factor)),
      
      params         = list(whichnitro        = whichnitro,
                            bgfactor          = bg.inflate.factor,
                            count.nbinom.size = count.nbinom.size,
                            verbose           = FALSE))
  }

RenderAllNitroExample <- 
  function(bg.inflate.factor = 1, count.nbinom.size = 20) {
    for (mynitro in c("NDEA", "NDMA", "NPIP", "NPYR")) {
      RenderNitroExample(mynitro, 
                         bg.inflate.factor = bg.inflate.factor,
                         count.nbinom.size = count.nbinom.size)
    }
  }



if (FALSE) {
  mSigAct:::RenderNitroExample("NDEA", bg.inflate.factor = 8)
}
