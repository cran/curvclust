
setMethod(
          f          = "EstepFMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD){
            
            Hcoef = 2^-(CCD@jklevels$scale*.Object@eta)
            
            lambda.Theta   = .Object@varE/.Object@Gamma2.Theta
            lambda.Nu      = .Object@varE/.Object@Gamma2.Nu
            
            Theta       = t(sapply(CCD@wavecoef,FUN = function(x){(x$D-.Object@Beta)/(1+lambda.Theta/Hcoef)}))
            Nu          = t(t(sapply(CCD@wavecoef,FUN = function(x){(x$C-.Object@Alpha)/(1+lambda.Nu)})))                
            .Object@Blup  = list(Theta=Theta,Nu=Nu)
            
            vardens.C           = .Object@Gamma2.Nu+.Object@varE
            vardens.D           = .Object@Gamma2.Theta*Hcoef+.Object@varE
            .Object@varBlup.Nu    = .Object@varE/(1+lambda.Nu)
            .Object@varBlup.Theta = .Object@varE/(1+lambda.Theta/Hcoef)
            
            logdens   = sapply(CCD@wavecoef, FUN = function(x){
              sum(dnorm( c(x$C,x$D),c(.Object@Alpha,.Object@Beta) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE))
              })            
            
            .Object@loglik = sum(logdens)
            return(.Object)
          }
          )
