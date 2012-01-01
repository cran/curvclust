\name{getwr.Y-methods}
\docType{methods}
\alias{getwr.Y-methods}
\alias{getwr.Y,CClustData-method}
\alias{getwr.Y}
\title{Inverse Wavelet Transform}
\description{Inverse Wavelet Transform for signal}
\usage{getwr.Y(.Object)}
\arguments{
\item{.Object}{Object of class CClustData}
}
\examples{
#library(curvclust)

## creation of a simulated dataset
#mu1          = c(rep(0,20),rep(1,20),rep(0,150),rep(-0.5,10),rep(1,60),rep(0.5,50),
# rep(0,200),rep(-1,80),rep(0,190),rep(-0.5,84),rep(0,80),rep(0.5,30),rep(0,50))
#mu2          = c(rep(0,20),rep(0.5,40),rep(0,130),rep(-0.5,10),rep(1,45),rep(0.5,65),
# rep(0,180),rep(-0.5,121),rep(0,86),rep(-0.5,170),rep(0,50),rep(1,57),rep(0,50))
#mu0          = cbind(mu1,mu2)
#n            = 50
#prop         = c(0.6,0.4)
#eta          = 2
#filter.number = 1

#SNR     = 5
#lambda  = 1
#simdata = simcurv(mu0,prop,SNR,lambda,n,eta,filter.number)
#simdata$CCD : simulated data
#simdata$CCR : parameters of the simulation

#CCD = simdata$CCD

## retrieve functional data
#Yfun =getwr.Y(CCD)
# plot functional data for individual 1
# plot(Yfun[[1]],type="l")

# dimension reduction is not performed to minimise reconstruction errors
# these plots illustrate that the reconstruction of the functional data 
# when coef. have been selected based on the unionset does not correspond to 
# traditional thresholding/reconstruction procedures

## dimension reduction
# CCDred = getUnionCoef(CCD)
# Yfun.red = getwr.Y(CCDred)
# plot(Yfun.red[[1]],type="l")


}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}
\keyword{methods}

