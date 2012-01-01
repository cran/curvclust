
setMethod(
          f          = "EstepFCM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD){
            
            logdens   = sapply(CCD@wavecoef, FUN = function(x){
              apply( dnorm( c(x$C,x$D),sapply(1:CCO@nbclust,FUN = function(ell){c(.Object@Alpha[[ell]],.Object@Beta[[ell]])}) , sd=sqrt(.Object@varE),log=TRUE),2,sum)                 
            })
            
            tau1           = log(.Object@prop)+logdens #n.mixt
            tau.max        = apply(tau1,2,max)#n
            tau1           = exp(t(tau1)-tau.max) #n*n.mixt            
            .Object@loglik = sum(log( apply(tau1,1,sum)) + tau.max)         
            .Object@Tau    = t(apply(tau1,1,FUN=function(x) x/sum(x)))
                                        
            return(.Object)
          }
          )
