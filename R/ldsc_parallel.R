#' 
#' @description 
#'
#' @
#'
ldsc_parallel <- function() {
  # Defaulting to PSOCK cluster as it should work on both Linux and Windows and from my experience it's faster when not copying large amounts of data
  cl <- parallel::makeCluster(int, type="PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  #Util-functions have to be explicitly passed to the analysis function in PSOCK cluster
  utilfuncs <- list()
  utilfuncs[[".LOG"]] <- .LOG
  
  foreach::foreach (i=1:length(filenames), .export=c(".munge_main"), .packages=c("stringr")) %dopar% {
    
  }
}

