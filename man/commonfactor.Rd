% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commonfactor.R
\name{commonfactor}
\alias{commonfactor}
\title{Run common factor model on genetic covariance S and sampling covariance matrix V}
\usage{
commonfactor(covstruc, estimation = "DWLS")
}
\arguments{
\item{covstruc}{Output from the multivariate LDSC/HDL/metaCCA function}

\item{estimation}{Options are either Diagonally Weighted Least Squares ("DWLS"; the default) or Maximum Likelihood ("ML")}
}
\value{
Common factor structure list
}
\description{
Function to run a common factor model based on output from multivariate LDSC/HDL/metaCCA
}
