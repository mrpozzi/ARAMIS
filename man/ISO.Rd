\name{ISO-class}
%\Rdversion{1.1}
\docType{class}
\alias{ISO-class}
\alias{[[,ISO-method}
\alias{$,ISO-method}
\alias{show,ISO-method}
\alias{weights,ISO-method}
\alias{mean,ISO-method}
\alias{var,ISO-method}
\alias{plot,ISO-method}
\alias{names,ISO-method}
\alias{attach,ISO-method}
\alias{sample,ISO-method}
\alias{summary,ISO-method}
\alias{ISO}
\alias{mean}
\alias{var}
\alias{plot}
\alias{sample}

\title{Class "ISO"}
\description{
An Importance Sampling Object
}
\section{Objects from the Class}{
The suggested way to construct objects of class \code{ISO} is to use
the constructor functions \code{AIS(N, niter, p, target, initialize, mixture, verbose = FALSE, seed = NULL)} and \code{AMIS(N, niter, p, target, initialize, mixture, verbose = FALSE, seed = NULL)}
}

\section{Slots}{
  \describe{
    \item{\code{IS}:}{A \code{"matrix"}  the Importance Sample in the first \code{p} columns and the Importance Weights in the last column.}
    \item{\code{Prop}:}{A \code{"numeric"} with the mixture proportion for the maximum perplexity mixture.}
    \item{\code{Mean}:}{A \code{"numeric"}
      containing the mean of the importance sample.}
    \item{\code{Var}:}{A \code{"matrix"}
       containing the Variance Covariance Matrix of the importance sample.}
    \item{\code{Perp}:}{A \code{"numeric"} with the maximum perplexity.}
    \item{\code{ESS}:}{A \code{"numeric"}, the Effective Sample Size.} 
\item{\code{envir}:}{A \code{"environment"}, it contains functions with suffix \code{r} and \code{p} with the "usual" \code{R} meaning, for simulating from the target and for computing the value of the CDF. Their Argument lists are:
\describe{
    \item{\code{rTarg(n)}}{}
    \item{\code{pTarg(q,lower.tail=TRUE,log.p=FALSE,N=1000)}}{}
}
} 
    \item{\code{seed}:}{A \code{"integer"}, the initial seed for the pseudo-random number generator.} 
    \item{\code{call}:}{The \code{"call"} to the function.}
    \item{\code{args}:}{A \code{"list"} with the arguments the function has been called with.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{AIS(N,niter,p,target,proposal=mvtComp(df=3),initialize,mixture,verbose=FALSE,tol=0.001,seed=sample(.Random.seed,1))
}}{A constructor
      function.  See also \code{\link{AIS}}}
    \item{\code{AMIS(N,niter,p,target,proposal=mvtComp(df=3),initialize=uniInit(),mixture,verbose=FALSE,parallel = c("no", "multicore", "snow"),nCores=-1,cl=NULL,tol=0.001,seed=sample(.Random.seed,1),...)}}{A constructor
      function.  See also \code{\link{AMIS}}}

\item{\code{plot(x,whichOne=1L,prj=1L:2L,N=100,xlim=c(-30,30),ylim=c(-30,30),main,...)}}{\code{signature(x = "ISO")}: plotting method.
\describe{
\item{\code{whichOne}}{Work in Progress}
\item{\code{prj}}{For \code{p}>2 a vector of length 2 specifying the 2D marginal to be plotted.}
\item{\code{N}}{Size of the sample to plot.}
}
}

\item{\code{mean(x, fun=NULL,...)}}{\code{signature(x = "ISO")}:  Returns the mean of \code{fun(x)}.
}

\item{\code{var(x, y=NULL,na.rm=FALSE,use)}}{\code{signature(x = "ISO")}: Returns the Variance-Covariance of \code{y(x)}.
}

\item{\code{weights}}{\code{signature(object = "ISO")}: Returns the IS weights.}

%\item{\code{attach}}{\code{signature(what = "ISO")}: \code{attach}es the %\code{envier} of the \code{ISO} object. It contains the following functions (substitute "Targ" with the name of the input \code{target}):
%\describe{
%\item{\code{pTarg(q,lower.tail=TRUE,log.p=FALSE,N=1000)}}{To compute the CDF of the Target}
%\item{\code{rTarg(IS)}}{To sample from the target distribution}
%}
%}

\item{\code{sample}}{\code{signature(x="ISO",size)}: To sample from the target distribution.}
\item{\code{names}}{\code{signature(x = "ISO")}: alias for \code{slotNames}.}


\item{\code{[[}}{\code{signature(x = "ISO")}: subsetting of the object.}
\item{\code{$}}{\code{signature(x = "ISO")}: subletting of the object.}

\item{\code{show}}{\code{signature(object = "ISO")}: the show method.}
  }
}

\author{
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
}
\seealso{
  The package vignette \code{demoARAMIS}.
}
\examples{
showClass("ISO")
}
\keyword{classes}
