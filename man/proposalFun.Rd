\name{proposal}
\alias{proposal}
\alias{mvtComp}
\alias{mvnormComp}

\title{
Proposal Distributions
}
\description{
Functions returning a component of a mixture for a proposal distribution. These functions provide an example of implementation for two interesting targets. They use a function closure approach in order to accept a general set of arguments and to store in their environment the variables specified by the user. The returned function has to be in the form

\code{
dproposal(xx,mu=rep(0,ncol(xx)),Sig=diag(1,ncol(xx),ncol(xx)))
}
\code{
rproposal(n,mu=0,Sig=1)
}


}

\arguments{
  \item{xx}{
A matrix with \code{p} columns
}
  \item{n}{
The size of the sample to simulate.
}
  \item{mu}{
\code{p} dimensional mean vector.
}
  \item{Sig}{
\code{p}X\code{p} Variance Covariance matrix.
}
}

\section{Available Proposals}{
\describe{
    \item{\code{mvtComp(df=3)}:}{Mixture of multivariate t distributions:
    \describe{
    \item{\code{df}}{ number of degrees of freedom.}
    }}
\item{\code{mvnormComp()}:}{Gaussian Mixture.}
}}

\value{
Value of the density for each observation.
A sample of size \code{n}.
}
\author{
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
}


\seealso{
See also \code{\link{target}}.
}
\examples{
mvtComp(df=3)$d(t(c(1,1)))
mvnormComp()$r(1)
}
