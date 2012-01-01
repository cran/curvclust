setMethod(
          f          = "getCM",
          signature  = "CClustData",
          definition =
          function(.Object,CCO) {
               
            CCR   = new("CClustRes",CCO,.Object)            
            n          = length(.Object@wavecoef)
            M          = length(.Object@wavecoef[[1]]$D) + 1
            Wd         = t(sapply(.Object@wavecoef,FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc         = t(t(sapply(.Object@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            
            CCR@Beta   = apply(Wd,2,mean)
            CCR@Alpha  = mean(Wc)        
            varE.D   = apply(Wd,2,var)*(n-1)
            varE.C   = var(Wc)*(n-1)
            CCR@varE   = c( (sum(varE.D) + varE.C )/ (n*M) )
            logdens    = sapply(.Object@wavecoef, FUN = function(x){
              sum(dnorm( c(x$C,x$D),c(CCR@Alpha,CCR@Beta), sd=sqrt(CCR@varE),log=TRUE))
            })
            CCR@loglik = sum(logdens)         
            return(CCR) 
          }
          )
