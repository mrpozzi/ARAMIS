\name{AIS-AMIS}
\alias{AIS}
\alias{AMIS}
\title{
Adaptive Importance Sampling
Adaptive Multiple Importance Sampling
}
\description{
Functions implementing Adaptive IS and constructors for object \code{ISO}.
}
\usage{
AIS(N,niter,p,target,proposal=mvtComp(df=3),initialize,mixture,verbose=FALSE,tol=0.001,seed=NULL)

AMIS(N,niter,p,target,proposal=mvtComp(df=3),initialize=uniInit(),mixture,verbose=FALSE,parallel = c("no", "multicore", "snow"),nCores=-1,cl=NULL,tol=0.001,seed=NULL,...)
}

\arguments{
  \item{N}{
Vector of length 3 with the sample size for the initialization and for the first and second Sampling phase, or a scalar if the three above are the same.
}
  \item{niter}{
Vector of length 2 with the number of iterations for the first and the second Sampling phase, or a scalar if the two above are the same.
}
  \item{p}{
Dimension of the sample space.
}
  \item{target}{
Target distribution. Takes as an argument a matrix with \code{p} columns and returns a value of the log-density for each row (observation). See \code{\link{target}}.
}
\item{proposal}{
A list with two components: \code{d} and \code{r}, to compute the Proposal Importance distribution and to simulate from it. Takes as an argument a matrix with \code{p} columns (or the size of the sample), a \code{p} dimensional vector \code{mu} and a \code{p}X\code{p} matrix \code{Sig}  and returns a value of the density for each observation (or the simulated sample). See \code{\link{proposal}}. Default is a 3 degrees of freedom multivariate t distribution.
}
  \item{initialize}{
See \code{\link{initialize}}.
Function for the initialization. Takes as an argument:
\describe{
\item{N0}{Sample Size of the Initialization phase}
\item{p}{Dimension of the sample space.}
\item{target}{As described above}
\item{verbose}{}
}
and returns a list with components:
\describe{
\item{w}{Vector with initial Importance weights.}
\item{xx}{Vector with initial Sample.}
\item{var}{Initial variance.}
\item{targ}{Initial Value of the target distribution (log-scale)}
\item{prop}{Initial Value of the proposal distribution (log-scale)}
}
Default is a random initialization.
}
  \item{mixture}{
See \code{\link{mixture}}.
Function for the parametric clustering. Takes as an argument a matrix with \code{p} columns and returns:
\describe{
\item{alpha}{Vector of Cluster Proportions}
\item{muHat}{Matrix with cluster means by rows. Dimensions \code{G}X\code{p}.}
\item{SigmaHat}{Array of Variance-Covariance matrix. Dimensions \code{p}X\code{p}X\code{G}.}
\item{G}{Number of components}
\item{cluster}{Labels vector.}
}
}
  \item{verbose}{
Default \code{FALSE}. Helps to spend the time while the function is running ;).
}
\item{parallel}{\code{"multicore"} doesn't work on Windows systems}
\item{nCores}{Number of cores to be used.}
\item{cl}{A \code{cluster} object to be used if by the prallel implementation.}
  \item{seed}{
initial seed for the (pseudo-)random number generator.
}
  \item{tol}{ a tolerance for the error}
  \item{...}{further arguments to be passed to subroutines}
}
\value{
An Object of class \code{ISO}, see \code{\link{ISO}}.
}
\references{
Jean-Marie Cornuet, Jean-Michel Marin, Antonietta Mira and Christian Robert (2012), Adaptive Multiple Importance Sampling, Scandinavian Journal of Statistics
}
\author{
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{initialize}}, \code{\link{mixture}} and \code{\link{target}}. See \code{demoARAMIS} vignette for more information.  
}
\examples{
ais <- AIS(N=50,niter=c(5,10),p=2,target= targetBanana(),initialize=uniInit(),mixture= mclustMix())
amis <- AMIS(N=c(10,50,100),niter=10,p=2,target=targetBanana(),initialize= amisInit(maxit=5), mixture=mclustMix(),verbose=TRUE)
}