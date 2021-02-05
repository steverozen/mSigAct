source("data-raw/errors/ID.attribution.error/code/debug-ID-attribution.R")

MAP.out <- DebugIDAttribution()

plus.minus.grid <- function(spectrum, sigs) {
  lapply(rownames(spectrum), 
        function(nn) { c(spectrum[nn, ] / sigs[nn, ]) }
        
        )
  
  
}
