\name{mixture}
\alias{mixture}
\alias{mclustMix}
\alias{fitMclust}
\alias{mixmodMix}


\title{
Model Based Clustering Routines
}
\description{
Functions for the fitting of the mixture model. They represent an example for the user to write his own. The returned function must be in the form

\code{
mixture(xx)
}
}


\arguments{
\item{xx}{
 A matrix with \code{p} columns 
}
}
\details{
\code{mixmodMix} is still work in progress...
}
\value{
Returns a list with components:
\item{alpha}{Vector of Cluster Proportions}
\item{muHat}{Matrix with cluster means by rows. Dimensions \code{G}X\code{p}.}
\item{SigmaHat}{Array of Variance-Covariance matrix. Dimensions \code{p}X\code{p}X\code{G}.}
\item{G}{Number of components.}
\item{cluster}{Labels vector.}
}

\author{
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
}

\examples{
# mclustMix(G = 2:9)(cbind(rnorm(100),rnorm(100)))
# mixmodMix(G = 1)(cbind(rnorm(100),rnorm(100)))
}
