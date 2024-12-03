#' Use 'snowfall' R package to parallelly estimate a genetic covariance structure using multivariate LD score regression
#' @description Function to run LD score regression (https://github.com/bulik/ldsc) to compute the genetic covariance between a series of traits based on genome wide summary statistics obtained from GWAS.
#' @references Consortium, S. W. G. of the P. G. et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat. Genet. 47, 291–295 (2015).
#' 
#' @param traits A vector of strings which point to LDSC munged files for trait you want to include in a SEM model
#' @param sample.prev A vector of sample prevalence for dichotomous traits and NA for continuous traits
#' @param population.prev A vector of population prevalence for dichotomous traits and NA for continuous traits
#' @param ld String which contains the path to the folder in which the LD scores used in the analysis are located.
#' @param wld String which contains the path to the folder in which the LD score weights used in the analysis are located
#' @param trait.names A character vector specifying how the traits should be named in the genetic covariance matrix (i.e., S). These variable names can subsequently be used in later steps for model specification. If no value is provided, the function will automatically name the variables using the generic from of V1-VX.
#' @param sep_weights Logical which indicates whether the weights are different form the LD scores used for the regression, defaults to FALSE.
#' @param chr defaults to 22
#' @param n.blocks defaults to 200
#' @param ldsc.log What to name the .log file if you want to override. Default to name the .log file based on file names used as input.
#' @param stand If you want to output the standard errors (SE) of the ld-score regression in the order they are listed in the genetic covariance matrix (i.e., S). Default to TRUE.
#' @param select default to FALSE, option: EVEN, ODD, c(...)
#' @param chisq.max default to NA
#' @param save.output default to FALSE, which indicates whether the output of ldsc function is saved to a .RData file named as ldsc.log.
#' @param output.path default to NULL (the current working directory). If you want to specify other output.path, please make output.path string ending with '/'.
#' @param cores default to 1
#'
#' @return The function returns a list with 5 named entries
#' @return S  estimated genetic variance/covariance matrix
#' @return SE variance matrix of the estimated parameter in S
#' @return I  matrix containing the "cross trait intercepts", or the error covariance between traits induced by overlap, in terms of subjects, between the GWAS on which the analyses are based
#' @return N  a vector containing the sample size for the genetic variances (i.e., heritability) and the geometric mean of sample sizes (i.e., sqrt(N1*N2)) between two samples for the genetic covariance
#' @return m  number of SNPs used to compute the LD scores
#' 
#' @export
ldsc_snowfall <- function(traits, sample.prev, population.prev, ld, wld,
                          trait.names = NULL, sep_weights = FALSE, chr = 22, n.blocks = 200,
                          ldsc.log = NULL, stand = TRUE, select = FALSE, chisq.max = NA,
                          save.output = FALSE, output.path = NULL, cores = 1) {
  
  # check traits
  if(length(traits) == 1) {warning("Our version of ldsc requires 2 or more traits. Please include an additional trait.")}
  
  # check trait.names
  if(is.null(trait.names)) {trait.names <- paste0("V", 1:length(traits))}
  if(!(is.null(trait.names))) {
    check_names <- stringr::str_detect(trait.names, "-")
    if(any(check_names))
      warning("Your trait names specified include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits using the trait.names argument.")
  }
  
  # check ldsc.log
  if(is.null(ldsc.log)) {
    logtraits <- gsub(".*/", "", traits)
    log2 <- paste(logtraits, collapse = "_")
    if (object.size(log2) > 200) {log2 <- substr(log2, 1, 100)}
    log.file <- file(paste0(log2, "_ldsc.log"), open="wt")
    output.filename <- paste0(log2, "_ldsc_output.RData")
  } else {
    log.file <- file(paste0(ldsc.log, "_ldsc.log"), open="wt")
    output.filename <- paste0(ldsc.log, "_ldsc_output.RData")
  }
  
  # check select
  if(select == "ODD" | select == "EVEN") {
    odd <- seq(1, chr, 2)
    even <- seq(2, chr, 2)
  }
  
  # check cores
  if (cores > length(traits)) {
    .LOG("Number of requested cores(", cores, ") greater than the number of munged files (", length(traits), "). Deferring to the lowest number", file=log.file)
    cores <- length(traits)
  }
  if (cores > parallel::detectCores()-1) {cores <- parallel::detectCores()-1}
  
  begin.time <- Sys.time()
  .LOG("Multivariate ld-score regression of ", length(traits), " traits ", "(", paste(traits, collapse = " "), ")", " began at: ", begin.time, file=log.file)
  
  # Dimensions of S and V matrix
  n.traits <- length(traits)
  n.V <- n.traits * (n.traits + 1) / 2
  
  if(n.traits > 18) {
    n.blocks <- (((n.traits + 1) * (n.traits + 2)) / 2) + 1
    n.blocks <- min(n.blocks, 1500)
    .LOG("     ", file=log.file, print=FALSE)
    .LOG("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to ", n.blocks, file=log.file)
    .LOG("This reflects the need to estimate V using at least one more block then they are nonredundant elements in the genetic covariance matrix that includes individual SNPs.", file=log.file)
    .LOG("If the n.blocks is > 1000 you should carefully inspect output for any strange results, such as extremely significant Q_SNP estimates.", file=log.file)
    .LOG("     ", file=log.file, print=FALSE)
    if(n.blocks > 1000) {
      warning("The number of blocks needed to estimate V is > 1000, which may result in sampling dependencies across the blocks used to estimate standard errors and can bias results.")
    }
  }
  
  #### Vector and Matrix Storage
  Liab.S <- rep(1, n.traits)
  cov <- matrix(NA, nrow=n.traits, ncol=n.traits)
  SE <- matrix(NA, nrow=n.traits, ncol=n.traits)
  I <- matrix(NA, nrow=n.traits, ncol=n.traits)
  N.vec <- matrix(NA, nrow=1, ncol=n.V)
  # V.hold <- matrix(NA, nrow=n.blocks, ncol=n.V)
  # V.hold <- bigmemory::big.matrix(nrow=n.blocks, ncol=n.V, type="double", init=NULL, separated=FALSE, backingfile=NULL, backingpath='/Volumes/Elements', descriptorfile=NULL, binarydescriptor=FALSE, shared=TRUE)
  
  #########  READ LD SCORES:
  .LOG('\n', "Reading in LD scores:", file=log.file)
  
  if(select == FALSE) {
    x <- do.call("rbind", lapply(1:chr, function(i) {
      suppressMessages(readr::read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
    .LOG(1:chr, ".l2.ldscore.gz", " have been read.", file=log.file)
  }
  
  if(select == "ODD") {
    x <- do.call("rbind", lapply(odd, function(i) {
      suppressMessages(readr::read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
    .LOG(odd, ".l2.ldscore.gz", " have been read.", file=log.file)
  }
  
  if(select == "EVEN") {
    x <- do.call("rbind", lapply(even, function(i) {
      suppressMessages(readr::read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
    .LOG(even, ".l2.ldscore.gz", " have been read.", file=log.file)
  }
  
  if(is.numeric(select)) {
    x <- do.call("rbind", lapply(select, function(i) {
      suppressMessages(readr::read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
    .LOG(select, ".l2.ldscore.gz", " have been read.", file=log.file)
  }
  
  x$CM <- NULL
  x$MAF <- NULL
  
  ######### READ weights:
  .LOG('\n', "Reading in weights of LD scores:", file=log.file)
  
  if(sep_weights) {
    if(select == FALSE) {
      w <- do.call("rbind", lapply(1:chr, function(i) {
        suppressMessages(readr::read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
      .LOG(1:chr, ".l2.ldscore.gz", " have been read.", file=log.file)
    }
    
    if(select == "EVEN") {
      w <- do.call("rbind", lapply(even, function(i) {
        suppressMessages(readr::read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
      .LOG(even, ".l2.ldscore.gz", " have been read.", file=log.file)
    }
    
    if(select == "ODD") {
      w <- do.call("rbind", lapply(odd, function(i) {
        suppressMessages(readr::read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
      .LOG(odd, ".l2.ldscore.gz", " have been read.", file=log.file)
    }
    
    if(is.numeric(select)) {
      w <- do.call("rbind", lapply(select, function(i) {
        suppressMessages(readr::read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
      .LOG(odd, ".l2.ldscore.gz", " have been read.", file=log.file)
    }
  } else {
    w <- x
    .LOG("Weights of LD scores and LD scores are identical.", file=log.file)
  }
  
  # rename last column of w, (i.e., L2 -> wLD)
  colnames(w)[ncol(w)] <- "wLD"
  w$CM <- NULL
  w$MAF <- NULL
  
  ###### READ M:
  .LOG('\n', "Reading in # of SNPs:", file=log.file)
  
  if(select == FALSE) {
    m <- do.call("rbind", lapply(1:chr, function(i) {
      suppressMessages(readr::read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
    .LOG(1:chr, ".l2.M_5_50", " have been read.", file=log.file)
  }
  
  if(select == "EVEN") {
    m <- do.call("rbind", lapply(even, function(i) {
      suppressMessages(readr::read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
    .LOG(even, ".l2.M_5_50", " have been read.", file=log.file)
  }
  
  if(select == "ODD") {
    m <- do.call("rbind", lapply(odd, function(i) {
      suppressMessages(readr::read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
    .LOG(odd, ".l2.M_5_50", " have been read.", file=log.file)
  }
  
  if(is.numeric(select)) {
    m <- do.call("rbind", lapply(select, function(i) {
      suppressMessages(readr::read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
    .LOG(select, ".l2.M_5_50", " have been read.", file=log.file)
  }
  
  m <- M.tot <- sum(m)
  
  
  # count the number of read summary statistics file
  s <- 0
  ## readSumstat function to read all summary statistics files and merge with LD-score files
  readSumstat <- function(chi1) {
    y1 <- suppressMessages(na.omit(readr::read_delim(chi1, delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)))
    # .LOG("Read in summary statistics [", s <<- s + 1, "/", n.traits, "] from: ", chi1, file=log.file)
    snowfall::sfCat(paste0("Read in summary statistics [", s <<- s + 1, "/", n.traits, "] from: ", chi1, '\n'))
    
    ## Merge with LD-score files
    merged <- merge(y1[, c("SNP", "N", "Z", "A1")], w[, c("SNP", "wLD")], by = "SNP", sort = FALSE)
    merged <- merge(merged, x, by = "SNP", sort = FALSE)
    merged <- merged[with(merged, order(CHR, BP)), ]
    # .LOG("Out of ", nrow(y1), " SNPs, ", nrow(merged), " remain after merging with LD-score files", file=log.file)
    snowfall::sfCat(paste0("Out of ", nrow(y1), " SNPs, ", nrow(merged), " remain after merging with LD-score files", '\n'))
    
    ## REMOVE SNPs with excess chisq.max:
    if(is.na(chisq.max)) {chisq.max <- max(0.001 * max(merged$N), 80)}
    rm <- (merged$Z^2 > chisq.max)
    merged <- merged[!rm, ]
    # .LOG("Removing ", sum(rm), " SNPs with Chi^2 > ", chisq.max, "; ", nrow(merged), " remain.", file=log.file)
    snowfall::sfCat(paste0("Removing ", sum(rm), " SNPs with Chi^2 > ", chisq.max, "; ", nrow(merged), " remain.", '\n'))
    
    # merged[SNP, N, Z, A1, wLD, CHR, BP, L2]
    return(merged)
  }
  
  .LOG('\n', "Parallelly reading in all summary statistics files:", file=log.file)
  ## Default to PSOCK cluster as it should work on both Linux and Windows and it's faster when not copying large amounts of data
  snowfall::sfInit(parallel = TRUE, cpus = cores, slaveOutfile = ldsc.log)
  # snowfall::sfLibrary(MASS)
  # snowfall::sfExport('')
  all_y <- snowfall::sfLapply(x=traits, fun=readSumstat)
  snowfall::sfStop()
  .LOG("All summary statistics files are parallelly read in.\n", file=log.file)
  return(all_y)
  
  
  # count the total number of runs in a double loop
  s <- 0
  
  for(j in 1:n.traits) {
    
    if(is.null(trait.names)) {chi1 <- traits[j]} else {chi1 <- trait.names[j]}
    
    y1 <- all_y[[j]]
    y1$chi1 <- y1$Z^2
    
  }
  
  ### Scale S to liability:
  ratio <- tcrossprod(sqrt(Liab.S))
  S <- cov * ratio
  
  ## Scale V to N per study (assume m constant)
  ## use crossprod instead of tcrossprod because N.vec is a one-row matrix
  # v.out <- cov(V.hold) / crossprod(N.vec * (sqrt(n.blocks) / m))
  ## calculate the ratio of the re-scaled and original S matrices
  # scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
  ## re-scale the sampling correlation matrix by the appropriate diagonals
  # V <- v.out * tcrossprod(scaleO)
  
  ## name traits according to trait.names argument and use general format of V1-VX if no trait.names provided
  colnames(S) <- if (is.null(trait.names)) paste0("V", 1:ncol(S)) else trait.names
  
  if(mean(Liab.S) != 1) {
    
    .LOG(c("     ", "     "), file=log.file, print=FALSE)
    .LOG("Liability Scale Results", file=log.file)
    
    for(j in 1:n.traits) {
      if(is.null(trait.names)) {chi1 <- traits[j]} else {chi1 <- trait.names[j]}
      for(k in j:length(traits)) {
        if(j == k) {
          .LOG("     ", file=log.file, print=FALSE)
          .LOG("Liability scale results for: ", chi1, file=log.file)
          .LOG("Total Liability Scale h2: ", round(S[j,k], 4), " (", round(SE[j,k], 4), ")", file=log.file)
        }
        if(j != k) {
          if(is.null(trait.names)) {chi2 <- traits[k]} else {chi2 <- trait.names[k]}
          .LOG("Total Liability Scale Genetic Covariance between ", chi1, " and ", chi2, ": ", round(S[k,j], 4), " (", round(SE[k,j], 4), ")", file=log.file)
          .LOG("     ", file=log.file, print=FALSE)
        }
      }
    }
  }
  
  end.time <- Sys.time()
  total.time <- difftime(time1=end.time, time2=begin.time, units="sec")
  mins <- floor(floor(total.time) / 60)
  secs <- floor(total.time - mins * 60)
  
  .LOG("     ", file=log.file, print=FALSE)
  .LOG('\n', "LDSC finished running at ", end.time, file=log.file)
  .LOG("Running LDSC for all files took ", mins, " minutes and ", secs, " seconds", file=log.file)
  .LOG("     ", file=log.file, print=FALSE)
  
  if(all(diag(S) > 0)) {
    
    if(stand) {
      ## calculate standardized results to print genetic correlations to log and screen
      ratio <- tcrossprod(1 / sqrt(diag(S)))
      S_Stand <- S * ratio
      
      ## Scale diagonal of matrix V to liability:
      V.diag <- gdata::lowerTriangle(SE, diag = TRUE)
      # calculate the ratio of the re-scaled and original S matrices
      scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
      # Make sure that if ratio in NaN (division by zero) we put the zero back in (In fact, not possible because of 'all(diag(S) > 0)')
      scaleO[is.nan(scaleO)] <- 0
      # re-scale the sampling correlation matrix by the appropriate diagonals
      V_Stand.diag <- V.diag * scaleO
      
      ## read unstandardized SEs from V.diag
      r <- nrow(S)
      SE_Stand <- matrix(0, r, r)
      SE_Stand[lower.tri(SE_Stand, diag = TRUE)] <- V_Stand.diag
      SE_Stand[upper.tri(SE_Stand, diag = FALSE)] <- t(SE_Stand)[upper.tri(SE_Stand, diag = FALSE)]
      
      .LOG(c("     ", "     "), file=log.file, print=FALSE)
      .LOG('\n', "Genetic Correlation Results", file=log.file)
      
      for(j in 1:n.traits) {
        if(is.null(trait.names)) {chi1 <- traits[j]} else {chi1 <- trait.names[j]}
        for(k in j:n.traits) {
          if(j != k) {
            if(is.null(trait.names)) {chi2 <- traits[k]} else {chi2 <- trait.names[k]}
            .LOG("Genetic Correlation between ", chi1, " and ", chi2, ": ", round(S_Stand[k, j], 4), " (", round(SE_Stand[k, j], 4), ")", file=log.file)
            .LOG("     ", file=log.file, print=FALSE)
          }
        }
      }
      
      flush(log.file)
      close(log.file)
      
      ldsc.output <- list(S=S, SE=SE, I=I, N=N.vec, m=m, S_Stand=S_Stand, SE_Stand=SE_Stand)
      if(save.output) {save(ldsc.output, file = paste0(output.path, output.filename))}
      return(list(S=S, SE=SE, I=I, N=N.vec, m=m, S_Stand=S_Stand, SE_Stand=SE_Stand))
    } else {
      ldsc.output <- list(S=S, SE=SE, I=I, N=N.vec, m=m)
      if(save.output) {save(ldsc.output, file = paste0(output.path, output.filename))}
      return(list(S=S, SE=SE, I=I, N=N.vec, m=m))
    }
    
  } else {
    .LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.", file=log.file)
    .LOG("Genetic correlation results could not be computed due to negative heritability estimates.", file=log.file)
    
    if(stand) {
      warning("Your genetic covariance matrix includes traits estimated to have a negative heritability. Please manually calculate standardized results. (e.g. You can use function Matrix::as.matrix((nearPD(x=ldsc.output$S, corr = TRUE))$mat))")
    }
    
    ldsc.output <- list(S=S, SE=SE, I=I, N=N.vec, m=m)
    if(save.output) {save(ldsc.output, file = paste0(output.path, output.filename))}
    return(ldsc.output)
  }
}
