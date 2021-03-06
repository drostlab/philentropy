% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{k_divergence}
\alias{k_divergence}
\title{K-Divergence (lowlevel function)}
\usage{
k_divergence(P, Q, testNA, unit)
}
\arguments{
\item{P}{a numeric vector storing the first distribution.}

\item{Q}{a numeric vector storing the second distribution.}

\item{testNA}{a logical value indicating whether or not distributions shall be checked for \code{NA} values.}

\item{unit}{type of \code{log} function. Option are 
\itemize{
\item \code{unit = "log"}
\item \code{unit = "log2"}
\item \code{unit = "log10"}   
}}
}
\description{
The lowlevel function for computing the k_divergence distance.
}
\examples{
k_divergence(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
}
\author{
Hajk-Georg Drost
}
