#' Clean and munge files to enable LD score regression
#' @description Function to process GWAS summary statistics files and prepare them for LD score regression
#' 
#' @param files A list of file names, files must be located in the working directory, or a path must be provided.
#' @param hm3 A file (UNZIPPED) of HAPMAP3 SNPs with some basic cleaning applied it is supplied and created by the original LD score regression developers and available here: https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
#' @param trait.names A vector of trait names which will be used as names for the munged files
#' @param N default to NULL
#' @param info.filter Numeric value which is used as a lower bound for imputation quality (INFO)
#' @param maf.filter Numeric value used as a lower bound for minor allele frequency
#' @param log.name What to name the .log file if you want to override. Default to name file based on file names used as input.
#' @param column.names column.names is a list (SNP, A1, A2, effect, P, MAF, Z) to specify other name
#' @param parallel default to FALSE
#' @param cores default to NULL
#' @param overwrite default to TRUE
#' @param output.path default to NULL (the current working directory). If you want to specify other output.path, please make output.path string ending with '/'.
#'
#' @return The function writes files of the ".sumstats.gz" format, which can be used to estimate SNP heritability and genetic covariance using the ldsc() function. The function will also output a .log file that should be examined to ensure that column names are being interpret correctly.
#' @export
munge <- function(files, hm3, trait.names = NULL, N = NULL, info.filter = .9, maf.filter = 0.01, log.name = NULL,
                  column.names = list(), parallel = FALSE, cores = NULL, overwrite = TRUE, output.path = NULL) {
  
  # check files
  if (is.list(files)) {
    warning(paste0("DeprecationWarning: In future versions a list of filenames will no longer be accepted.\n",
                  "                    Please change files to a vector to ensure future compatibility."))
    files_ <- c()
    for (i in 1:length(files)) {files_ <- c(files_, files[[i]])}
    files <- files_
  }
  filenames <- as.vector(files)
  
  # check N
  if (is.null(N)) {N <- rep(NA, length(files))}
  
  # Sanity checks start
  .check_file_exists(files)
  .check_file_exists(hm3)
  .check_equal_length(files, trait.names)
  .check_equal_length(files, N)
  .check_range(N, min=0, max=Inf, allowNA=TRUE)
  .check_range(info.filter)
  .check_range(maf.filter)
  if (any(!(names(column.names) %in% c("SNP", "A1", "A2", "effect", "INFO", "P", "N",  "MAF", "Z")))) {
    stop(paste0("Names in column.names not recognized. Please use the following keys:\n        ",
                paste(c("SNP", "A1", "A2", "effect", "INFO", "P", "N", "MAF", "Z"), collapse=", ")))
  }
  .check_boolean(parallel)
  if (!is.null(cores)) .check_range(cores, min=0, max=Inf)
  .check_boolean(overwrite)
  # Sanity checks finished

  # check log.name
  if(is.null(log.name)) {
    log2 <- paste(trait.names, collapse="_")
    if(nchar(log2) > 200) {log2 <- substr(log2, 1, 100)}
    log.file <- file(paste0(log2, "_munge.log"), open="wt")
  } else {
    log.file <- file(paste0(log.name, "_munge.log"), open = "wt")
  }
  
  # check overwrite
  existing_traits <- c()
  if (overwrite) {
    existing_files <- c()
    for (trait.name in trait.names) {
      if (file.exists(paste0(trait.name, ".sumstats.gz"))) {
        existing_traits <- c(existing_traits, trait.name)
        existing_files <- c(existing_files, paste0(trait.name, ".sumstats.gz"))
      }
    }
    if (length(existing_files) > 0)
      .LOG("File(s) ", paste0(existing_files, collapse = ", "), " already exist and will be overwritten", file=log.file)
  } else {
    trait.names <- trait.names[-which(trait.names == existing_traits)]
  }
  
  # Read HAPMAP3 SNPs reference file
  .LOG("Reading in reference file", file=log.file)
  ref <- data.table::fread(hm3, header=T, data.table=F)
  
  begin.time <- Sys.time()
  .LOG("The munging of ", length(trait.names), " summary statistics started at ", begin.time, file=log.file)
  
  # check parallel
  if (parallel & (length(files) == 1)) {
    .LOG("Parallel munging requested for a single file.\nParallel munging only has benefits for munging multiple files.\nParallel disabled", file=log.file)
    parallel <- FALSE
  }
  if (!parallel) {
    .LOG("Reading summary statistics for ", paste(files, collapse=" "), ". Please note that this step usually takes a few minutes due to the size of summary statistic files.", file=log.file)
    ## note that fread is not used here due to formatting differences across summary statistic files
    files <- lapply(files, read.table, header=T, quote="\"", fill=T, na.string=c(".", NA, "NA", ""))
    .LOG("All summary statistic files are loaded into R!", file=log.file)
    for(i in 1:length(files)) {
      .munge_main(i, NULL, files[[i]], filenames[i], trait.names[i], N[i], ref, hm3, info.filter, maf.filter, column.names, overwrite, output.path, log.file)
    }
  } else {
    if(is.null(cores)) {int <- parallel::detectCores() - 1}
    else {int <- cores}
    if (int > length(filenames)) {
      .LOG("Number of requested cores(", int, ") greater than the number of files (", length(filenames), "). Deferring to the lowest number", file=log.file)
      int <- length(filenames)
    }
    
    # Defaulting to PSOCK cluster and it's faster when not copying large amounts of data
    cl <- parallel::makeCluster(int, type="PSOCK")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    #Util-functions have to be explicitly passed to the analysis function in PSOCK cluster
    utilfuncs <- list()
    utilfuncs[[".get_renamed_colnames"]] <- .get_renamed_colnames
    utilfuncs[[".LOG"]] <- .LOG
    utilfuncs[["gzip"]] <- R.utils::gzip
    .LOG("As parallel munging was requested, logs of each sumstats file will be saved separately", file=log.file)
    foreach::foreach (i=1:length(filenames), .export=c(".munge_main"), .packages=c("stringr")) %dopar% {
      .munge_main(i, utilfuncs, NULL, filenames[i], trait.names[i], N[i], ref, hm3, info.filter, maf.filter, column.names, overwrite, output.path, NULL)
    }
  }
  
  end.time <- Sys.time()
  total.time <- difftime(time1=end.time, time2=begin.time, units="sec")
  mins <- floor(floor(total.time) / 60)
  secs <- total.time - mins*60
  
  .LOG("     ", file=log.file)
  .LOG("Munging was completed at ", end.time, file = log.file)
  .LOG("The munging of all files took ", mins, " minutes and ", secs, " seconds", file=log.file)
  .LOG("Please check the .log file(s) to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files", file=log.file)
  
  flush(log.file)
  close(log.file)
}