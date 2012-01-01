\name{CClustRes-class}
\Rdversion{1.1}
\docType{class}
\alias{CClustRes-class}
\alias{initialize,CClustRes-method}
\alias{summary,CClustRes-method}
\alias{[,CClustRes-method}
\alias{[<-,CClustRes-method}

\title{Class "CClustRes"}
\description{Class containing the Results for Curve Clustering}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("CClustRes", CCO, CCD)}.
\describe{
	\item{\code{CCO}:}{Object of class \code{CClustO} (options)}
	\item{\code{CCD}:}{Object of class \code{CClustD} (data)}
	 } 
}
\section{Slots}{
  \describe{
    \item{\code{prop}:}{\code{"numeric"} or \code{"NULL"}, proportions of the mixture model}
    \item{\code{Beta}:}{\code{"list"} or \code{"numeric"}, wavelet coefficients of the fixed effects}
    \item{\code{Alpha}:}{\code{"list"} or \code{"numeric"}, scaling coefficients of the fixed effects}
    \item{\code{Blup}:}{\code{"list"} or \code{"NULL"}, Best Linear Predictors of the functional random effects}
    \item{\code{varBlup.Nu}:}{\code{"list"} or \code{"numeric"} or \code{"NULL"}, variance of the Best Linear Predictors of the scaling coefficients of the random effects}
    \item{\code{varBlup.Theta}:}{\code{"list"} or \code{"numeric"} or \code{"NULL"}, variance of the Best Linear Predictors of the wavelet coefficients of the random effects}
    \item{\code{varE}:}{\code{"numeric"}, error variance}
    \item{\code{Gamma2.Theta}:}{\code{"list"} or \code{"numeric"} or \code{"NULL"},variance of the scaling coefficients of the random effects}
    \item{\code{Gamma2.Nu}:}{ \code{"list"} or \code{"numeric"} or \code{"NULL"},variance of the wavelet coefficients of the random effects}
    \item{\code{Tau}:}{\code{"matrix"} or \code{"NULL"}, matrix of posterior probabilities of cluster membership}
    \item{\code{eta}:}{\code{"numeric"} or \code{"NULL"}, estimated regularity of the functional random effect}
    \item{\code{loglik}:}{\code{"numeric"}, loglikelihood of the model}
  }
}
\section{Methods}{
  \describe{
	\item{initialize}{}
	\item{[}{}
	\item{summary}{summary of CCR object}
	\item{getBIC}{Calculates the BIC from estimated parameters}
	\item{getICL}{Calculates the ICL from estimated parameters}
	\item{getwr.mu}{Returns the inverse wavelet transform for mean curves mu}
	\item{getwr.Blup}{Returns the inverse wavelet transform for functional random effects predictions (functional BLUPs)}
	\item{getwr.varU}{Returns the estimated variance/covariance matrix of the functional random effects}
	\item{getlambda.U}{Calculates the estimated variance ratio between error variances and random effect variances}
	\item{getSNR.mu}{Calculates the Signal to Noise Ratio for functional fixed effects}
	 }
}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}
\examples{
showClass("CClustRes")

library(curvclust)
showClass("CClustData")

# creation of a simulated dataset
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


#CCD = simdata$CCD

## dimension reduction
#CCDred = getUnionCoef(CCD)

## options setting
#CCO            = new("CClustO")
#CCO["nbclust"] = 2

## computing Functional Clustering Mixed Model
#CCR         = getFCMM(CCDred,CCO)
#summary(CCR)

}
\keyword{classes}
