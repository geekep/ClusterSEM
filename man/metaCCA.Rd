% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaCCA.R
\name{metaCCA}
\alias{metaCCA}
\title{Estimate a genetic covariance matrix using using canonical correlation analysis based on the work by Cichonska, A. et al., Bioinformatics, 2016}
\usage{
metaCCA(traits, sample.prev = NA, population.prev = NA, trait.names = NULL)
}
\arguments{
\item{traits}{A vector of strings which point to munged files for trait you want to include in a SEM model. The HDL function works with standard munged files.}

\item{sample.prev}{A vector of sample prevalence for dichotomous traits and NA for continuous traits}

\item{population.prev}{A vector of population prevalence for dichotomous traits and NA for continuous traits}

\item{trait.names}{A character vector specifying how the traits should be named in the genetic covariance matrix (i.e., S). These variable names can subsequently be used in later steps for model specification. If no value is provided, the function will automatically name the variables using the generic from of V1-VX.}
}
\value{
The function returns a list with 3 named entries

S	estimated genetic variance/covariance matrix

SE variance matrix of the estimated parameter in S

I matrix containing the "cross trait intercepts", or the error covariance between traits induced by overlap, in terms of subjects, between the GWAS on which the analyses are based
}
\description{
Function to run metaCCA (https://github.com/aalto-ics-kepaco/metaCCA-matlab) to compute the genetic covariance between a series of traits based on genome-wide summary statistics.
}
\references{
Cichonska, A. et al. metaCCA: summary statistics-based multivariate meta-analysis of genome-wide association studies using canonical correlation analysis. Bioinformatics 32, 1981–1989 (2016).
}
