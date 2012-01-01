
setMethod(
          f          = "EstepFCMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD){
            
            
            Hcoef = 2^-(CCD@jklevels$scale*.Object@eta)
            
            if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){
              
              lambda.Theta = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE/.Object@Gamma2.Theta[[ell]]})
              lambda.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE/.Object@Gamma2.Nu[[ell]]})
              
              .Object@Blup  = sapply(1:CCO@nbclust, FUN=function(ell){
                Theta = t(sapply(CCD@wavecoef,FUN = function(x){
                  (x$D-.Object@Beta[[ell]])/(1+lambda.Theta[[ell]]/Hcoef)
                }))
                Nu    = t(t(sapply(CCD@wavecoef,FUN = function(x){
                  (x$C-.Object@Alpha[[ell]])/(1+lambda.Nu[[ell]])
                })))
                return(list(Theta=Theta,Nu=Nu))
              },simplify=FALSE,USE.NAMES=TRUE)
              
              vardens.C           = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Nu[[ell]]+.Object@varE})                
              vardens.D           = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Theta[[ell]]*Hcoef+.Object@varE})
              .Object@varBlup.Nu    = lapply(1:CCO@nbclust,FUN=function(ell){.Object@varE/(1+lambda.Nu[[ell]])})              
              .Object@varBlup.Theta = lapply(1:CCO@nbclust,FUN=function(ell){.Object@varE/(1+lambda.Theta[[ell]]/Hcoef)})
              if (CCO@Gamma2.structure=="group"){.Object@varBlup.Theta = lapply(.Object@varBlup.Theta,FUN=function(y){mean(y)})}        	              
              	             
              
            }
            
            
            if ((CCO@Gamma2.structure=="scale.location") | (CCO@Gamma2.structure=="constant")) {
              
              lambda.Theta = .Object@varE/.Object@Gamma2.Theta
              lambda.Nu    = .Object@varE/.Object@Gamma2.Nu
              
              .Object@Blup  = sapply(1:CCO@nbclust, FUN=function(ell){
                Theta = t(sapply(CCD@wavecoef,FUN = function(x){
                  (x$D-.Object@Beta[[ell]])/(1+lambda.Theta/Hcoef)
                }))
                Nu    = t(t(sapply(CCD@wavecoef,FUN = function(x){
                  (x$C-.Object@Alpha[[ell]])/(1+lambda.Nu)
                })))
                return(list(Theta=Theta,Nu=Nu))
              },simplify=FALSE,USE.NAMES=TRUE)
              
              vardens.C           = .Object@Gamma2.Nu+.Object@varE
              vardens.D           = .Object@Gamma2.Theta*Hcoef+.Object@varE      
              .Object@varBlup.Nu    = .Object@varE/(1+lambda.Nu)
              .Object@varBlup.Theta = .Object@varE/(1+lambda.Theta/Hcoef)
              if (CCO@Gamma2.structure=="constant"){
              	.Object@varBlup.Theta = mean(.Object@varBlup.Theta)
              }
              
            }
            
            logdens   = sapply(CCD@wavecoef, FUN = function(x){
              apply( dnorm( c(x$C,x$D),sapply(1:CCO@nbclust,FUN = function(ell){c(.Object@Alpha[[ell]],.Object@Beta[[ell]])}) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE),2,sum)
            })
                        
            tau1           = log(.Object@prop)+logdens 
            tau.max        = apply(tau1,2,max)
            tau1           = exp(t(tau1)-tau.max)
            .Object@loglik = sum(log( apply(tau1,1,sum)) + tau.max)         
            .Object@Tau    = t(apply(tau1,1,FUN=function(x) x/sum(x)))
            
            return(.Object)
          }
          )
