#'  Re-sampling from the genetic covariance matrix and variance covariance matrix of the parameter estimates in S (i.e., LDSC/HDL/metaCCA)
#' @description Function for re-sampling from the genetic covariance matrix and variance covariance matrix of the parameter estimates in S
#' 
#' @param suffStat Output from the multivariate LDSC/HDL/metaCCA function
#' @param seed default to NULL
#' @param rep default to 250
#'
#' @return Re-sampled genetic covariance matrix, # of matrix = rep
#' @export
ClusterSEM.suffStat <- function(suffStat,
                          seed = NULL,
                          rep = 250) {
  set.seed(seed = seed)
  if (is.null(suffStat$S)) {
    baz <-
      MASS::mvrnorm(rep, suffStat$S_Full[lower.tri(suffStat$S_Full, diag = T)], Sigma = suffStat$V_Full)
    suffStat$baz <- baz
  } else {
    baz <-
      MASS::mvrnorm(rep, suffStat$S[lower.tri(suffStat$S, diag = T)], Sigma = suffStat$V)
    suffStat$baz <- baz
  }
  
  suffStat$rep <- rep
  return(suffStat)
}

#' Independence test for genetic covariance structure
#' @description Function to independence test for genetic covariance structure
#'
#' @param suffStat Output from the multivariate LDSC/HDL/metaCCA function
#' @param x An item from 1:nrow(suffStat$S)
#' @param y An item from 1:nrow(suffStat$S)
#' @param S A set from 1:nrow(suffStat$S)
#'
#' @return p value indicating independence between x and y
#' @export
ClusterSEM.indepTest <- function(x, y, S, suffStat) {
  if (is.null(suffStat$S)) {
    sig <- matrix(NA, nrow(suffStat$S_Full), ncol(suffStat$S_Full))
  } else {
    sig <- matrix(NA, nrow(suffStat$S), ncol(suffStat$S))
  }
  
  pcor <- matrix(NA, suffStat$rep, 1)
  
  for (i in 1:suffStat$rep) {
    sig[lower.tri(sig, diag = T)] <- suffStat$baz[i,]
    sig <- makeSymm(sig)
    pcor[i,] <-
      -corpcor::cor2pcor(sig[c(x, y, S), c(x, y, S)])[1, 2]
  }
  
  p_val <-
    pchisq((mean(pcor) / stats::sd(pcor)) ^ 2, 1, lower.tail = F)
  return(p_val)
}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
