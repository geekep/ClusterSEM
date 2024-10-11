#' Automate writing lavaan model syntax based on clustering output
#' @description Function to automate writing lavaan model syntax based on clustering output
#'
#' @param S_LD Genetic covariance matrix
#' @param fix_resid default to TRUE
#' @param common Whether to specify a common factor model
#'
#' @return
#' @export
write.model <-
  function(clusters,
           S_LD,
           fix_resid = TRUE,
           common = FALSE) {
    Model <- ""
    if (common == TRUE) {
      for (f in 1) {
        u <- 1
        Model1 <- ""
        for (i in 1:nrow(S_LD)) {
          if (u == 1) {
            linestart <- paste("F", f, "=~",  colnames(S_LD)[i], sep = "")
            u <- u + 1
            linemid <- ""
          } else {
            linemid <- paste(linemid, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      Model <- paste(Model, linestart, linemid, " \n ", sep = "")
    } else {
      for (f in 1:ncol(clusters)) {
        u <- 1
        Model1 <- ""
        for (i in 1:nrow(Loadings)) {
          if (u == 1) {
            linestart <- paste("F", f, "=~",  colnames(S_LD)[i], sep = "")
            u <- u + 1
            linemid <- ""
          } else{
            linemid <- paste(linemid, " + ", colnames(S_LD)[i], sep = "")
          }
          
        }
        Model <- paste(Model, linestart, linemid, " \n ", sep = "")
        linestart <- ""
        linemid <- ""
      }
    }
    
    if (fix_resid == TRUE) {
      Model3 <- ""
      #create unique combination of letters for residual variance parameter labels
      n <-
        combn(letters, 4)[, sample(1:14000, ncol(S_LD), replace = FALSE)]
      for (i in 1:ncol(S_LD)) {
        if (grepl(colnames(S_LD)[i], Model) == TRUE) {
          linestart3a <-
            paste(colnames(S_LD)[i],
                  " ~~ ",
                  paste(n[, i], collapse = ""),
                  "*",
                  colnames(S_LD)[i],
                  sep = "")
          linestart3b <-
            paste(paste(n[, i], collapse = ""), " > .0001", sep = "")
          Model3 <-
            paste(Model3, linestart3a, " \n ", linestart3b, " \n ", sep = "")
        }
      }
      Model <- paste(Model, Model3, sep = "")
    }
    
    return(Model)
  }
