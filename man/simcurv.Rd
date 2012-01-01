\name{simcurv}
\alias{simcurv}
\title{simcurv}
\description{Simulation of clusters of curves}
\arguments{
  \item{mu}{Cluster average function normalized by their L2 norm (matrix
  of size M*L)}
  \item{prop}{Proportions of each cluster}
  \item{SNR}{Signal to Noise Ratio}
  \item{lambda}{Ratio between error variance and random effect variances}
  \item{n}{Number of curves (individuals) to simulate}
  \item{eta}{Random process regularity}
  \item{filter.number}{Filter of the wavelet transform for the DaubExPhase family (see wd imported from the wavethresh package}
  \item{Gamma2.structure}{Structure of the random effect variance (constant, scale.location)}

}
\value{
  \item{out}{List containing an object of class CClustData (simulated
    data) and an object of class CClustRes to store the true parameters
    of the simulation}
}
\examples{

## creation of a simulated dataset (only one group)
#mu0          = c(rep(0,20),rep(1,20),rep(0,150),rep(-0.5,10),rep(1,60),rep(0.5,50),
# rep(0,200),rep(-1,80),rep(0,190),rep(-0.5,84),rep(0,80),rep(0.5,30),rep(0,50))
#mu0          = matrix(mu0)
#n            = 50
#prop         = 1
#eta          = 2
#filter.number= 1

#SNR     = 5
#lambda  = 1
#simdata = simcurv(mu0,prop,SNR,lambda,n,eta,filter.number)

#simdata$CCD : simulated data
#simdata$CCR : parameters of the simulation

#CCD = simdata$CCD


## creation of a simulated dataset with 2 groups
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


}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}
