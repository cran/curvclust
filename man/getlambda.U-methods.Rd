\name{getlambda.U-methods}
\docType{methods}
\alias{getlambda.U-methods}
\alias{getlambda.U,CClustRes-method}
\alias{getlambda.U}
\title{get lambda.U}
\description{Computes the ratio between the estimated variance of the functional random effects and the estimated noise variance}
\usage{getlambda.U(.Object,CCO,CCD)}
\arguments{
\item{.Object}{Object of class CClustRes}
\item{CCO}{Object of class CClustO}
\item{CCD}{Object of class CClustData}
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

## dimension reduction
#CCDred = getUnionCoef(CCD)


## options setting
#CCO            = new("CClustO")
#CCO["nbclust"] = 1

## computing Functional Mixed Model
#CCR         = getFMM(CCDred,CCO)
#summary(CCR)
#getlambda.U(CCR,CCO)

## options setting
#CCO            = new("CClustO")
#CCO["nbclust"] = 2

## computing Functional Clustering Mixed Model
#CCR         = getFCMM(CCDred,CCO)
#summary(CCR)
#getlambda.U(CCR,CCO,CCDred)


}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{Madison Giacofci, Sophie Lambert-Lacroix, Guillemette Marot, Franck Picard}
\keyword{methods}

