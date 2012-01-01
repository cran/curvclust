\name{getFMM-methods}
\docType{methods}
\alias{getFMM-methods}
\alias{getFMM,CClustData-method}
\alias{getFMM}
\title{get Functional Mixed Model}
\description{Estimates Parameters for a Functional Mixed Model without clustering (one group)}
\usage{getFMM(.Object,CCO)}
\arguments{
\item{.Object}{Object of class CClustData}
\item{CCO}{Object of class CClustO}
}

\examples{
## creation of a simulated dataset (only one group)
#mu0         = c(rep(0,20),rep(1,20),rep(0,150),rep(-0.5,10),rep(1,60),rep(0.5,50),
# rep(0,200),rep(-1,80),rep(0,190),rep(-0.5,84),rep(0,80),rep(0.5,30),rep(0,50))
#mu0         = matrix(mu0)
#n           = 50
#prop        = 1
#eta         = 2
#filter.num  = 1

#SNR     = 5
#lambda  = 1
#simdata = simcurv(mu0,prop,SNR,lambda,n,eta,filter.num)

#simdata$CCD : simulated data
#simdata$CCR : parameters of the simulation

#CCD = simdata$CCD

## options setting
#CCO            = new("CClustO")
#CCO["nbclust"] = 1

## computing Functional Mixed Model
#CCR         = getFMM(CCD,CCO)
#summary(CCR)

}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}
\keyword{methods}

