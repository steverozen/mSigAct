test_small_MAPAssignActivity <- function() {

  # Create test data: 2x3 spectra matrix and 2x4 signatures matrix

  sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh38$SBS96[ , c("SBS1", "SBS2"), drop = FALSE]
  spect <- round(rowSums(sigs) * 1000)
  
    
  # Normalize signatures to sum to 1 (typical requirement for mutational signatures)
  sigs <- apply(sigs, 2, function(x) x / sum(x))
  
  # Create temporary output directory
  output.dir <- file.path(tempdir(), "test_small_MAP")
  
  browser()
  # Call MAPAssignActivity with minimal parameters
  result <- mSigAct:::MAPAssignActivity1(
    spect = spect,
    sigs  = sigs,
    seed = 456,
    drop.low.mut.samples = FALSE      # Don't drop samples for small test
  )
  browser()
  # Clean up
  unlink(output.dir, recursive = TRUE)
  
  # Return result for inspection
  return(result)
}

# Run the test
result <- test_small_MAPAssignActivity()
print("Test completed successfully!")
print(paste("Result structure:", paste(names(result), collapse = ", ")))