% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KL.R
\name{KL}
\alias{KL}
\title{Kullback-Leibler Divergence}
\usage{
KL(x, test.na = TRUE, unit = "log2", est.prob = NULL)
}
\arguments{
\item{x}{a numeric \code{data.frame} or \code{matrix} (storing probability vectors) or a numeric \code{data.frame} or \code{matrix} storing counts (if \code{est.prob = TRUE}). See \code{\link{distance}} for details.}

\item{test.na}{a boolean value indicating whether input vectors should be tested for NA values.}

\item{unit}{a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.}

\item{est.prob}{method to estimate probabilities from a count vector. Default: est.prob = NULL.}
}
\value{
The Kullback-Leibler divergence of probability vectors.
}
\description{
This function computes the Kullback-Leibler divergence of two probability
distributions P and Q.
}
\details{
\deqn{KL(P||Q) = \sum P(P) * log2(P(P) / P(Q)) = H(P,Q) - H(P)}

where H(P,Q) denotes the joint entropy of the probability
distributions P and Q and H(P) denotes the entropy of
probability distribution P. In case P = Q then KL(P,Q) = 0 and in case P !=
Q then KL(P,Q) > 0.

The KL divergence is a non-symmetric measure of the directed divergence
between two probability distributions P and Q. It only fulfills the
\emph{positivity} property of a \emph{distance metric}.

Because of the relation KL(P||Q) = H(P,Q) - H(P), the Kullback-Leibler
divergence of two probability distributions P and Q is also named
\emph{Cross Entropy} of two probability distributions P and Q.
}
\examples{

# Kulback-Leibler Divergence between P and Q
P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)
x <- rbind(P,Q)
KL(x)

# Kulback-Leibler Divergence between P and Q using different log bases
KL(x, unit = "log2") # Default
KL(x, unit = "log")
KL(x, unit = "log10")

# Kulback-Leibler Divergence between count vectors P.count and Q.count
P.count <- 1:10
Q.count <- 20:29
x.count <- rbind(P.count,Q.count)
KL(x, est.prob = "empirical")

# Example: Distance Matrix using KL-Distance

Prob <- rbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))

# compute the KL matrix of a given probability matrix
KLMatrix <- KL(Prob)

# plot a heatmap of the corresponding KL matrix
heatmap(KLMatrix)


}
\references{
Cover Thomas M. and Thomas Joy A. 2006. Elements of Information
Theory. \emph{John Wiley & Sons}.
}
\seealso{
\code{\link{H}}, \code{\link{CE}}, \code{\link{JSD}}, \code{\link{gJSD}}, \code{\link{distance}}
}
\author{
Hajk-Georg Drost
}
