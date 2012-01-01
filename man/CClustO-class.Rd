\name{CClustO-class}
\Rdversion{1.1}
\docType{class}
\alias{CClustO-class}
\alias{[,CClustO-method}
\alias{[<-,CClustO-method}

\title{Class "CClustO"}
\description{Class containing the Options for Curve Clustering}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("CClustO", ...)}}
\section{Slots}{
  \describe{
    \item{\code{nbclust}:}{\code{"numeric"}, number of clusters}
    \item{\code{Gamma2.structure}:}{\code{"character"} or \code{"NULL"}, structure of the variance of the random effect, (\code{"constant"}, \code{"group"}, \code{"scale.location"}, \code{"group.scale.location"}, \code{"none"})}
    \item{\code{burn}:}{\code{"numeric"}, number of burn iterations for initialization}
    \item{\code{eps}:}{\code{"numeric"}, precision (relative) of the EM algorithm}
    \item{\code{loglikplot}:}{\code{"logical"}, display of the log-likelihood function with the number of iterations of EM}
    \item{\code{init}:}{\code{"character"}, initialization method, \code{"SEM"} (stochastic EM) or \code{"rEM"} (random seed EM)}   
  }
}

\section{Methods}{
  \describe{
	\item{initialize}{}
	\item{[}{}
	 }
}  
\examples{
showClass("CClustO")

## options setting
#CCO            = new("CClustO")
# to set the number of clusters
#CCO["nbclust"] = 1
# to set the structure of the random effect
#CCO["Gamma2.structure"] = "constant"
# to consider clustering models without random effects
#CCO["Gamma2.structure"] = "none"
# to set the number of burn iterations of SEM and rEM
#CCO["burn"] = 100
# to set the (relative) precision of the EM algorithm stopping rule
#CCO["eps"] = 0.01
# to plot the likelihood along initialization + EM iterations
#CCO["loglikplot"] = TRUE
# to set the initialization method
#CCO["init"] = "SEM"
}
\references{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard, Wavelet-based clustering for mixed-effects functional models in high dimension}
\author{M. Giacofci, S. Lambert-Lacroix, G. Marot, F. Picard}

\keyword{classes}
