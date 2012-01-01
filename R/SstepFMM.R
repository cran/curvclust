setMethod(
          f          = "SstepFMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {            
            n              = length(CCD@wavecoef)            
            M              = length(CCD@wavecoef[[1]]$D) + 1
            Wd             = t(sapply(CCD@wavecoef,  FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc             = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef          = 2^-(CCD@jklevels[,1]*.Object@eta)
            Theta        = t(apply(.Object@Blup$Theta,1,FUN = function(x){rnorm(M-1,x,sqrt(Hcoef*.Object@varBlup.Theta))}))
            Nu           = t(t(apply(.Object@Blup$Nu,1,FUN = function(x){rnorm(1,x,sqrt(.Object@varBlup.Nu))})))
            .Object@Blup   = list(Theta = Theta,Nu = Nu)
            return(.Object)
          }
          )
