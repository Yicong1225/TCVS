% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/P.R
\docType{data}
\name{P}
\alias{P}
\title{Taxonomy Structure of OTUs in our example}
\format{
An object of class `matrix`.
}
\usage{
data(P_60)
}
\description{
The taxonomy structure of OTUs.
The first 2p rows of the matrix P correspond to the $l_1$ penalty 
applied to the CLR (Centered Log-Ratio) covariates followed by their knockoff copies. 
The remaining rows correspond to the tree structure associated with the CLR covariates and their knockoff copies, 
following the same order of covariates first and knockoff copies second.
}
\keyword{datasets}
