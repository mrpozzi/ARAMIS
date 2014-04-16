\name{ARAMIS-package}
\alias{ARAMIS-package}
\alias{ARAMIS}

\docType{package}
\title{
A R Adaptive Multiple Importance Sampling
}
\description{
ARAMIS is an R package that runs the AMIS algorithm in parallel. ARAMIS exploits the embarrassingly parallel structure of AMIS to improve the performance of the algorithm and reduce simulation time by orders of magnitudes (the reduction in computational time is directly proportional to the number of particles used in the underlying Importance Sampler). Following a Map-Reduce scheme the computations are sent to the independent computing units and then the results are collected and merged together. Each core takes care of an approximately equal set of particles evaluating the importance weights and merging them together to compute the new IS estimates at each iteration of the algorithm. The musketeers slogan "all for one, one for all" The potential of the R language is fully exploited by ARAMIS which allows the user to taylor the software to the model which results from his or her own research setting. Target and proposal kernel can be easily specified by the user. Some working examples contained in the manual explain how this can be efficiently and easily done.
As a consequence of the flexibility of the package even fairly complicated problems can be accommodated, e.g. sampling from a Copula distribution with a mixture of extreme value distribution as the proposal kernel. Very little ad-hoc code is needed: package \code{copula} provides a way to evaluate the target distribution and package \code{evd} offers useful tools for extreme value distributions which are easily plugged in the ARAMIS function.
}
\details{
\tabular{ll}{
Package: \tab AMIS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-01\cr
License: \tab GPL (>= 3)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
Antonietta Mira,
\email{antonietta.mira@usi.ch}

Maintainer:
Luca Pozzi,
\email{p.luc@stat.berkeley.edu}
}
\references{
Jean-Marie Cornuet, Jean-Michel Marin, Antonietta Mira and Christian Robert (2012), Adaptive Multiple Importance Sampling, Scandinavian Journal of Statistics
}
\keyword{ package }
\seealso{
See also the package vignette \code{demoARAMIS}.
}
