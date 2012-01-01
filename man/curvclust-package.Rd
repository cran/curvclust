\name{curvclust-package}
\alias{curvclust-package}
\alias{curvclust}
\docType{package}
\title{
Curve clustering
}
\description{
\code{curvclust} is a package dedicated to model-based curve clustering. Considered models include Functional Clustering Mixed Models (FCMM, ie functional clustering with the presence of functional random effects),
but also traditional functional clustering model (FCM, without functional random effects), and functional mixed models (FMM, functional random effects without clustering). Estimation is done by maximum likelihood using the EM algorithm, and two criteria are 
proposed to select the number of clusters, based on integrated likelihoods.
}
\details{
The main object classes are \code{"CClustData"},\code{"CClustO"}, \code{"CClustRes"} to store data, options and results respectively.

\code{getFCM}, \code{getFCMM} and \code{getFMM} are the most important methods for clustering. 

\code{simcurv} can be used to simulate signals.

Methods \code{getwr.Y}, \code{getwr.mu}, \code{getwr.Blup}, \code{getwr.varU} and \code{getwr.varBlup} perform 
inverse wavelet transforms to retrieve the functional observations, the functional fixed effects, the functional Blups, and the variance and conditional 
variance matrices of the random effect and Blup. 

See the corresponding help pages for more details.

}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}
\examples{
#showClass("CClustO")
#showClass("CClustData")
#showClass("CClustRes")

## creation of a simulated dataset
#mu1         = c(rep(0,20),rep(1,20),rep(0,150),rep(-0.5,10),rep(1,60),rep(0.5,50),
#rep(0,200),rep(-1,80),rep(0,190),rep(-0.5,84),rep(0,80),rep(0.5,30),rep(0,50))
#mu2         = c(rep(0,20),rep(0.5,40),rep(0,130),rep(-0.5,10),rep(1,45),rep(0.5,65),
#rep(0,180),rep(-0.5,121),rep(0,86),rep(-0.5,170),rep(0,50),rep(1,57),rep(0,50))

#simdata = simcurv(mu=cbind(mu1,mu2),prop=c(0.6,0.4),SNR=5,lambda=1,n=50,eta=2,
#	          filter.number=1)

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

#groups      = apply(CCR["Tau"],1,which.max)
}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}

\keyword{ package }
