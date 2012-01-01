\name{getFCM-methods}
\docType{methods}
\alias{getFCM-methods}
\alias{getFCM}
\alias{getFCM,CClustData-method}
\title{get Functional Clustering Model}
\description{Estimates Parameters for a Functional Clustering Model without random effect}
\usage{getFCM(.Object,CCO)}
\arguments{
\item{.Object}{Object of class CClustData}
\item{CCO}{Object of class CClustO}
}
\examples{
## creation of a simulated dataset
#mu1          = c(rep(0,20),rep(1,20),rep(0,150),rep(-0.5,10),rep(1,60),rep(0.5,50),
# rep(0,200),rep(-1,80),rep(0,190),rep(-0.5,84),rep(0,80),rep(0.5,30),rep(0,50))
#mu2          = c(rep(0,20),rep(0.5,40),rep(0,130),rep(-0.5,10),rep(1,45),rep(0.5,65),
# rep(0,180),rep(-0.5,121),rep(0,86),rep(-0.5,170),rep(0,50),rep(1,57),rep(0,50))
#mu0          = cbind(mu1,mu2)
#n           = 50
#prop        = c(0.6,0.4)
#eta         = 2
#filter.num  = 1

#SNR     = 5
#lambda  = 1
#simdata = simcurv(mu0,prop,SNR,lambda,n,eta,filter.num)


##simdata$CCD : simulated data
##simdata$CCR : parameters of the simulation

#CCD = simdata$CCD

## dimension reduction
#CCDred = getUnionCoef(CCD)

## options setting
#CCO         = new("CClustO")
#CCO["nbclust"] = 2
#CCO["Gamma2.structure"] = "none"

## computing Functional Clustering Model
#CCR         = getFCM(CCDred,CCO)
#summary(CCR)
}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}
\keyword{methods}

