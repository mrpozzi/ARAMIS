\name{target}
\alias{target}
\alias{targetBanana}
\alias{targetMVN}
\alias{targetMix}

\title{
Target Distributions
}

\description{
Functions returning a target distribution. These function provide an example of implementation for three interesting targets. They use a function closure approach in order to accept a general set of arguments and to store in their environment the variables specified by the user. The returned function must be in the form

\code{
target(x)
}


}


\arguments{
  \item{x}{
A \code{p}-dimensional observation
}
}

\section{Example of Target Distributions}{
\describe{
    \item{\code{targetBanana(sig=100,b=0.01,log=TRUE)}:}{Banana Shaped Target from Haario et al. (1999):
    \describe{
    \item{\code{sig}}{ Variance of x_2;}
    \item{\code{b}}{ "Banana-ness" parameter.
                     Moderate banana shape: b=0.01. 
                     Strong banana shape: b=0.03.}
    \item{log}{ logical: should the log-density be returned.}
    }}
\item{\code{targetMVN(Mu=rep(0,2),Sigma=diag(1,2,2))}:}{Multivariate Normal Target:
    \describe{
    \item{\code{Mu}}{ Vector of Means;}
    \item{\code{Sigma}}{ Covariance Matrix.}
    }}
\item{\code{targetMix(alpha=rep(.5,2),Mu=matrix(0,2,2),Sigma=array(cbind(diag(1,2,2),diag(1,2,2)),dim=c(2,2,2)))}:}{
Gaussian Mixture Target:
    \describe{
    \item{\code{alpha}}{ Vector of Mixture probabilities;}
    \item{\code{Mu}}{ Matrix of Means (by row);}
    \item{\code{Sigma}}{ Array of Covariance Matrices.}
    }
}
}}

\value{
The value of the log-density.
}
\references{
Haario, H., Saksman, E., and Tamminen, J. (1999). Adaptive proposal distribution for random walk Metropolis algorithm. Computational Statistics, 14:375-395.
}
\author{
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
}

\seealso{
See also \code{\link{proposal}}. For more details see the tutorial in \code{vignette("demoARAMIS")}.
}

\examples{
targetBanana()(c(1,2))
targetMVN()(c(1,2))
targetMix()(c(1,2))
}
