#' @keywords internal
NullReturnForMAPAssignActivity <- 
  function(signature.universe,
           target.spectrum = NULL,
           msg         = "Message not set", 
           time.used   = system.time(3),
           use.forward.search = FALSE) {
    
    if (is.null(target.spectrum)) {
      sample.name <- "No sample"
    } else {
      sample.name <- colnames(target.spectrum)
    }
    
    null.assignment <- matrix(rep(0, ncol(signature.universe)))
    colnames(null.assignment) <- sample.name
    rownames(null.assignment) <- colnames(signature.universe)
    
    null.spect       <- matrix(rep(0, nrow(signature.universe)))
    colnames(null.spect) <- sample.name
    rownames(null.spect)  <- rownames(signature.universe)
    if (!is.null(target.spectrum)) {
      null.spect <- CopyAttributes(from = target.spectrum, to = null.spect)
    }
    
    all.tested <- data.frame(
      sig.names        = "",
      p.for.sig.subset = NA,
      exp              = NA,
      loglh.of.exp     = NA,
      df               = NA,
      akaike.weights   = NA,
      LRT.p.value      = NA,
      LRT.q.value      = NA
    )
    
    alt.solutions <- # An 0-row data frame (no alt.solutions)
      structure(
        list(sig.names = character(0), 
             p.for.sig.subset = numeric(0),      
             exp = list(), 
             loglh.of.exp = numeric(0), 
             df = numeric(0),
             akaike.weights = logical(0),
             LRT.p.value = numeric(0),
             LRT.q.value = numeric(0),
             QP.exp = list(),
             QP.cosine = numeric(0)),
        row.names = integer(0), 
        class = "data.frame")
    
    default.distances <-
      structure(
        list(
          method = c("log.likelihood", "euclidean", "manhattan",  "cosine"),
          proposed.assignment = c(log.likelihood = NA,
                                  euclidean      = NA, 
                                  manhattan      = NA, 
                                  cosine         = NA),
          QP.assignment       = c(log.likelihood = NA,
                                  euclidean      = NA,
                                  manhattan      = NA,
                                  cosine         = NA)),
        class = c("tbl_df",  "tbl", "data.frame"), 
        row.names = c(NA, -4L))
    
    xx <- 
      list(proposed.assignment           = null.assignment, # Not sure if this should be a 0-row matrix but that is a scalar 0
           proposed.reconstruction       = null.spect,
           reconstruction.distances      = default.distances,
           time.for.assignment           = time.used,
           error.messages                = msg
      )
    
    if (use.forward.search) {
      return(xx)
    } else {
      return(c(xx, 
               list(all.tested      = all.tested,
                    alt.solutions   = alt.solutions)))
    }
  }