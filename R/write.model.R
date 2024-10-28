#' Automate writing lavaan model syntax based on clustering output
#' @description Function to automate writing lavaan model syntax based on clustering output
#' 
#' @param S_LD Genetic covariance matrix
#' @param clusters Output from various clustering algorithms
#' @param common Whether to specify a common factor model, default to FALSE
#' @param hierarchical Whether to specify a hierarchical exploratory factor model, default to FALSE
#' @param fix_resid default to TRUE
#'
#' @return lavaan model syntax satisfying clusters
#' @export
write.model <-
  function(S_LD,
           clusters=list(),
           common = FALSE,
           hierarchical = FALSE,
           fix_resid = TRUE) {
    
    #final model syntax satisfying clusters
    Model <- ""
    
    if (common == TRUE) {
      
      for (f in 1) {
        u <- 1
        Model1 <- ""
        for (i in 1:nrow(S_LD)) {
          if (u == 1) {
            line.start <- paste("F", f, "=~", colnames(S_LD)[i], sep = "")
            u <- u + 1
            line.mid <- ""
          } else {
            line.mid <- paste(line.mid, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      Model1 <- paste(Model1, line.start, line.mid, " \n ", sep = "")
      Model <- paste(Model, Model1, sep = "")
    } else if (hierarchical == TRUE) {
      
      Model2 <- ""
      
    } else {
      
      Model3 <- ""
      for (f in 1:length(clusters[["size"]])) {
        
        IDP.ID <- colnames(S_LD)[stringr::str_detect(clusters[[".Data"]], as.character(f), negate = FALSE)]
        
        u <- 1
        for (i in 1:length(IDP.ID)) {
          if (u == 1) {
            line.start <- paste("F", f, " =~ ", IDP.ID[i], sep = "")
            u <- u + 1
            line.mid <- ""
          } else{
            line.mid <- paste(line.mid, " + ", IDP.ID[i], sep = "")
          }
        }
        Model3 <- paste(Model3, line.start, line.mid, " \n ", sep = "")
        line.start <- ""
        line.mid <- ""
      }
      Model <- paste(Model, Model3, sep = "")
    }
    
    if (fix_resid == TRUE) {
      Model4 <- ""
      #create unique combination of letters for residual variance parameter labels
      n <- combn(letters, 4)[, sample(1:14000, ncol(S_LD), replace = FALSE)]
      for (i in 1:ncol(S_LD)) {
        if (grepl(colnames(S_LD)[i], Model) == TRUE) {
          residual.variance <-
            paste(colnames(S_LD)[i],
                  " ~~ ",
                  paste(n[, i], collapse = ""),
                  "*",
                  colnames(S_LD)[i],
                  sep = "")
          residual.variance.condition <-
            paste(paste(n[, i], collapse = ""), " > .0001", sep = "")
          Model4 <-
            paste(Model4, residual.variance, " \n ", residual.variance.condition, " \n ", sep = "")
        }
      }
      Model <- paste(Model, Model4, sep = "")
    }
    
    return(Model)
  }
