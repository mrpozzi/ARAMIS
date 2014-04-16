\name{initialize}
\alias{initialize}
\alias{amisInit}
\alias{varInit}
\alias{uniInit}

\title{
Initialization Routines
}
\description{
Functions returning a routine to initialize the IS algorithm. They use a function closure approach in order to accept a general set of arguments and to store in their environment the variables specified by the user.
The returned function must be in the form

\code{
initialize(N0,p,target,proposal,verbose=FALSE)
}
}

\arguments{
\item{N0}{Sample Size of the Initialization phase}
\item{p}{Dimension of the sample space.}
\item{target}{As described in \code{\link{target}}}
\item{proposal}{As described in \code{\link{proposal}}}
\item{verbose}{To receive updates on the initialization.}
}


\section{Initialization Procedure}{
\describe{
    \item{\code{amisInit(maxit=5000,maxVar=100,s=sqrt(maxVar))}}{Initialization procedure described in Cornuet et al. (2012):
    \describe{
    \item{\code{maxit}}{ Maximum number of iterations for the optimization routine;}
    \item{\code{maxVar}}{upper bound of the variance.}
    \item{s}{scale parameter.}
}
    }
\item{\code{varInit(Var)}:}{
    \describe{
    \item{\code{Var}}{Searches over a range of variances for the one that maximizes the ESS of a sample from the proposal.}
    }
}
\item{\code{uniInit()}}{ Random initialization.}
}
}

\value{
\describe{
Returns a list with components:
\item{w}{Vector with initial Importance weights.}
\item{xx}{Vector with initial Sample.}
\item{var}{Initial variance.}
\item{targ}{Initial Value of the target distribution (log-scale)}
\item{prop}{Initial Value of the proposal distribution (log-scale)}
}
}

\references{
Jean-Marie Cornuet, Jean-Michel Marin, Antonietta Mira and Christian Robert (2012), Adaptive Multiple Importance Sampling, Scandinavian Journal of Statistics
}

\author{
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
}

\examples{
N0 <- 5
p <- 2
proposal <- mvtComp(df=3)
target <- targetBanana()
amisInit()(N0,p,target,proposal)
varInit(Var=c(10,100))(N0,p,target,proposal)
uniInit()(N0,p,target,proposal)
}

