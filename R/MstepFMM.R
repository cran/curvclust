setMethod(
          f          = "MstepFMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {
            
            n              = length(CCD@wavecoef)
            M              = length(CCD@wavecoef[[1]]$D) + 1
            Wd             = t(sapply(CCD@wavecoef,FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc             = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef          = 2^-(CCD@jklevels[,1]*.Object@eta)
            
            .Object@Beta  = apply(Wd-.Object@Blup$Theta,2,mean)
            .Object@Alpha = mean(Wc-.Object@Blup$Nu)
            
            lambda.Theta = .Object@varE/.Object@Gamma2.Theta
            lambda.Nu    = .Object@varE/.Object@Gamma2.Nu
            
              
            varE.D     = (t(Wd-.Object@Blup$Theta) -.Object@Beta)^2 + .Object@varE/(1+lambda.Theta/Hcoef)
            varE.C     = sum( (Wc-.Object@Blup$Nu-.Object@Alpha)^2 + .Object@varE/(1+lambda.Nu))
            .Object@varE = c( (sum(varE.D)+ varE.C) / (n*M) )
            
            .Object@Gamma2.Nu = sum( .Object@Blup$Nu^2 + .Object@varE /(1+lambda.Nu) )/n
            Z                 = ( t(.Object@Blup$Theta^2) + .Object@varE /(1+lambda.Theta/Hcoef) ) / Hcoef
            
            if ((CCO@Gamma2.structure=="constant") | (CCO@Gamma2.structure=="group")) {
              .Object@Gamma2.Theta = sum(Z) / (n*(M-1))
            } else{
              .Object@Gamma2.Theta = apply(Z,1,mean)
            }
            .Object@eta=getetaFMM(.Object,CCD,CCO)
            return(.Object) 
          }
          )

setGeneric( name = "getllkFMM" ,def = function(.Object,CCD,CCO,eta){standardGeneric("getllkFMM")})
setMethod(
          f          = "getllkFMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCD,CCO,eta) {
                 
            Hcoef        = 2^-(CCD@jklevels$scale*eta)            
            lambda.Theta = .Object@varE/.Object@Gamma2.Theta
            lambda.Nu    = .Object@varE/.Object@Gamma2.Nu            
            Theta      = t(sapply(CCD@wavecoef,FUN = function(x){(x$D-.Object@Beta)/(1+lambda.Theta/Hcoef)}))
            Nu         = t(t(sapply(CCD@wavecoef,FUN = function(x){(x$C-.Object@Alpha)/(1+lambda.Nu)})))                
            vardens.C  = .Object@Gamma2.Nu+.Object@varE
            vardens.D  = .Object@Gamma2.Theta*Hcoef+.Object@varE            
            logdens   = sapply(CCD@wavecoef, FUN = function(x){
              sum(dnorm( c(x$C,x$D),c(.Object@Alpha,.Object@Beta) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE))
              })                        
            llk = sum(logdens)       
            return(llk)
          })


setGeneric( name = "getetaFMM" ,def = function(.Object,CCD,CCO){standardGeneric("getetaFMM")})
setMethod(
          f          = "getetaFMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCD,CCO) {
            
            x1  = 1+1e-4        
            x3  = 6
            r   = (sqrt(5)-1)/2
            x2  = x1+(1-r)*(x3-x1)
            x4  = x1+r*(x3-x1)
            Jx2 = -getllkFMM(.Object,CCD,CCO,x2)
            Jx4 = -getllkFMM(.Object,CCD,CCO,x4)

            while (abs(x4-x2)>10^-4){
              if (Jx2<Jx4){
                x3        = x4
                x4        = x2
                Jx4       = Jx2
                x2        = x1+(1-r)*(x3-x1)
                Jx2       = -getllkFMM(.Object,CCD,CCO,x2)

              } else {
                if (Jx2>Jx4){
                  x1         = x2
                  x2         = x4
                  Jx2        = Jx4
                  x4         = x1+r*(x3-x1)
                  Jx4       = -getllkFMM(.Object,CCD,CCO,x4)
                } else {
                  x1        = x2
                  x3        = x4
                  x2        = x1+(1-r)*(x3-x1)
                  x4        = x1+r*(x3-x1)        
                  Jx2       = -getllkFMM(.Object,CCD,CCO,x2)
                  Jx4       = -getllkFMM(.Object,CCD,CCO,x4)
                }
                
              }
            }
            etamax=max(round((x2+x4)/2,4),1+1e-4)
            etamax
          })
