\name{CClustData-class}
\Rdversion{1.1}
\docType{class}
\alias{CClustData-class}
\alias{initialize,CClustData-method}
\alias{print,CClustData-method}
\alias{show,CClustData-method}
\alias{summary,CClustData-method}
\alias{[,CClustData-method}
\alias{[<-,CClustData-method}


\title{Class "CClustData"}
\description{Class containing the Data for Curve Clustering}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("CClustData", Y, filter.number)}.
\describe{
    \item{\code{Y}:}{List containing the curves for each individual (\code{Y[[i]]} is the curve for individual \code{i}). List names are conserved in CClustData objects.}   
    \item{\code{filter.number}:}{Filter of the wavelet transform for the DaubExPhase family (see \code{wd} imported from the wavethresh package}
    }
}
\section{Slots}{
  \describe{
    \item{\code{wavecoef}:}{List of Wavelet coefficients (C: approximation coef, D: wavelet coef)}
    \item{\code{jklevels}:}{Scales (j) and positions (k) of the wavelet coefficients}
    \item{\code{seMAD}:}{Vector of MAD variance estimates for each curve}
    \item{\code{lengthsignal}:}{Length of the signal}
    \item{\code{filter.number}:}{Filter of the wavelet transform for the DaubExPhase family (see \code{wd} imported from the wavethresh package)}
    }
  }
\section{Methods}{
  \describe{
	\item{initialize}{}
	\item{[}{}
	\item{show}{}
	\item{summary}{Summary of CCD object}
	\item{getCM}{Computes a constant model}
	\item{getFCM}{Computes Functional Clustering Model}
	\item{getFMM}{Computes Functional Mixed Model}
	\item{getFCMM}{Computes Functional Clustering Mixed Model}
	\item{getUnionCoef}{Selects the Union of Wavelet Coefficients for Dimension Reduction}
	 }
}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}
\examples{
#library(curvclust)
#showClass("CClustData")

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
##simdata$CCD : simulated data
##simdata$CCR : parameters of the simulation

#CCD=simdata$CCD
#show(CCD)
#summary(CCD)

## dimension reduction
#CCDred = getUnionCoef(CCD)

#summary(CCDred)

## to create a CClustData object directly from data
#mu1          = c(rep(0,20),rep(1,20),rep(0,150),rep(-0.5,10),rep(1,60),rep(0.5,50),
# rep(0,200),rep(-1,80),rep(0,190),rep(-0.5,84),rep(0,80),rep(0.5,30),rep(0,50))
#mu2          = c(rep(0,20),rep(0.5,40),rep(0,130),rep(-0.5,10),rep(1,45),rep(0.5,65),
# rep(0,180),rep(-0.5,121),rep(0,86),rep(-0.5,170),rep(0,50),rep(1,57),rep(0,50))
#n    = 4
#M    = length(mu1)
##each entry of the list should be a curve for one individual
#Y    = list()
#Y[[1]] = mu1 + rnorm(M,0,0.1)
#Y[[2]] = mu1 + rnorm(M,0,0.1)
#Y[[3]] = mu2 + rnorm(M,0,0.1)
#Y[[4]] = mu2 + rnorm(M,0,0.1)
#matplot(t(Reduce("rbind",Y)),type="l")
#CCD    = new("CClustData",Y=Y,filter.number=1)

}
\keyword{classes}

