#' Load data from IEU GWAS summary dataset
#' @description Estimate the causal effects of body mass index || id:ieu-a-2 on coronary artery disease || id:ieu-a-7 risk using 79 uncorrelated SNP genetic variants
#' @param exposure default to "ieu-a-2"
#' @param outcome default to "ieu-a-7"
#'
#' @return harmonised data
#' 
#' @importFrom TwoSampleMR extract_instruments
#' @importFrom TwoSampleMR extract_outcome_data
#' @importFrom TwoSampleMR harmonise_data
#' @importFrom ieugwasr check_access_token
get_coef_from_gwas <- function(exposure = "ieu-a-2",
                               outcome  = "ieu-a-7") {
  
  a <-  extract_instruments(outcomes = exposure,
                           p1 = 5e-08,
                           clump = TRUE,
                           p2 = 5e-08,
                           r2 = 0.001,
                           kb = 10000,
                           access_token = check_access_token(),
                           force_server = TRUE)
  
  b <-  extract_outcome_data(snps = a$SNP,
                            outcomes = outcome,
                            proxies = TRUE,
                            rsq = 0.8,
                            align_alleles = 1,
                            palindromes = 1,
                            maf_threshold = 0.3,
                            access_token = check_access_token(),
                            splitsize = 10000,
                            proxy_splitsize = 500)
  
  dat <-  harmonise_data(exposure_dat = a,
                        outcome_dat = b,
                        action = 2)
  
  return(dat)
}


#' Load multi data from IEU GWAS summary dataset
#' @description Estimate the multivariate effects of HDL cholesterol || id:ieu-a-299 (#SNP=79), LDL cholesterol || id:ieu-a-300 (#SNP=68) and triglycerides || id:ieu-a-302 (#SNP=42) on coronary artery disease || id:ieu-a-7 risk using 143 uncorrelated SNP genetic variants
#' @param exposures default to c("ieu-a-299","ieu-a-300","ieu-a-302")
#' @param outcome default to "ieu-a-7"
#'
#' @return harmonised data
#' 
#' @importFrom TwoSampleMR mv_extract_exposures
#' @importFrom TwoSampleMR extract_outcome_data
#' @importFrom TwoSampleMR mv_harmonise_data
#' @importFrom ieugwasr check_access_token
mv_get_coef_from_gwas <- function(exposures = c("ieu-a-299",
                                                "ieu-a-300",
                                                "ieu-a-302"),
                                  outcome   = "ieu-a-7") {
  
  a <- mv_extract_exposures(id_exposure = exposures,
                            clump_r2 = 0.001,
                            clump_kb = 10000,
                            harmonise_strictness = 2,
                            access_token = check_access_token(),
                            find_proxies = TRUE,
                            force_server = TRUE,
                            pval_threshold = 5e-08,
                            pop = "EUR")
  
  b <-  extract_outcome_data(snps = a$SNP,
                            outcomes = outcome,
                            proxies = TRUE,
                            rsq = 0.8,
                            align_alleles = 1,
                            palindromes = 1,
                            maf_threshold = 0.3,
                            access_token = check_access_token(),
                            splitsize = 10000,
                            proxy_splitsize = 500)
  
  dat <-  mv_harmonise_data(exposure_dat = a,
                           outcome_dat = b,
                           harmonise_strictness = 2)
  
  return(dat)
}
