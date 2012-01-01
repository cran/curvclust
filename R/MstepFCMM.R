setMethod(
          f          = "MstepFCMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {

            ## avoid label switching, sort groups by levels
            tmp           = order(unlist(.Object@Alpha))
            .Object@Beta  = .Object@Beta[tmp]
            .Object@Alpha = .Object@Alpha[tmp]
            .Object@prop  = .Object@prop[tmp]
            .Object@Tau   = .Object@Tau[,tmp]
            .Object@Blup  = .Object@Blup[tmp]
            
            n             = length(CCD@wavecoef)
            M             = length(CCD@wavecoef[[1]]$D) + 1
            Wd            = t(sapply(CCD@wavecoef,FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc            = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef         = 2^-(CCD@jklevels[,1]*.Object@eta)

            clustfreq     = apply(.Object@Tau,2,sum)
            .Object@prop  = as.vector(clustfreq) / n 
            .Object@Beta  = lapply(1:CCO@nbclust,FUN = function(ell){
              apply( .Object@Tau[,ell]* (Wd-.Object@Blup[[ell]]$Theta),2,sum)/clustfreq[ell]
            })
            .Object@Alpha = lapply(1:CCO@nbclust,FUN = function(ell){
              sum(.Object@Tau[,ell]*(Wc-.Object@Blup[[ell]]$Nu))/clustfreq[ell]
            })

            if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){

              .Object@Gamma2.Theta = .Object@Gamma2.Theta[tmp]
              .Object@Gamma2.Nu    = .Object@Gamma2.Nu[tmp]
              
              lambda.Theta = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE /.Object@Gamma2.Theta[[ell]]})
              lambda.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){.Object@varE /.Object@Gamma2.Nu[[ell]]})
              
              varE.D =  sapply(1:CCO@nbclust,FUN = function(ell){
                RSS = apply(Wd-.Object@Blup[[ell]]$Theta,1,FUN = function(x){x-.Object@Beta[[ell]]})^2 + .Object@varE/(1+lambda.Theta[[ell]]/Hcoef) 
                sum(.Object@Tau[,ell]*t(RSS))
              },simplify=FALSE,USE.NAMES=TRUE)
              
              varE.C =  sapply(1:CCO@nbclust,FUN = function(ell){
                RSS = (Wc-.Object@Blup[[ell]]$Nu-.Object@Alpha[[ell]])^2 + .Object@varE/(1+lambda.Nu[[ell]])
                sum(.Object@Tau[,ell]*RSS)
              },simplify=FALSE,USE.NAMES=TRUE)
              
              .Object@varE = (Reduce("+",varE.D) + Reduce("+",varE.C)) / (n*M)
              
              .Object@Gamma2.Nu = sapply(1:CCO@nbclust,FUN = function(ell){
                Z = .Object@Blup[[ell]]$Nu^2 + .Object@varE /(1+lambda.Nu[[ell]])
                sum(.Object@Tau[,ell]*Z) /clustfreq[ell]
              },simplify=FALSE,USE.NAMES=TRUE)
              
              Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                Z = (t(.Object@Blup[[ell]]$Theta^2) + .Object@varE /(1+lambda.Theta[[ell]]/Hcoef) ) / Hcoef
                return(t(Z)*.Object@Tau[,ell])
              },simplify=FALSE,USE.NAMES=TRUE)

              if( (CCO@Gamma2.structure=="group.scale.location") ){
                
                .Object@Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                  apply(Gamma2.Theta[[ell]],2,sum)/clustfreq[ell]
                },simplify=FALSE,USE.NAMES=TRUE)
                
                
              } else if (CCO@Gamma2.structure=="group"){
                
                .Object@Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                  sum(Gamma2.Theta[[ell]])/(clustfreq[ell]*(M-1))
                },simplify=FALSE,USE.NAMES=TRUE)
                
              }
            }
            
            if ((CCO@Gamma2.structure=="scale.location") | (CCO@Gamma2.structure=="constant")) {
              
              lambda.Theta = .Object@varE/.Object@Gamma2.Theta
              lambda.Nu    = .Object@varE/.Object@Gamma2.Nu                             
              
              varE.D =  sapply(1:CCO@nbclust,FUN = function(ell){
                RSS = apply(Wd-.Object@Blup[[ell]]$Theta,1,FUN = function(x){x-.Object@Beta[[ell]]})^2 + .Object@varE/(1+lambda.Theta/Hcoef) 
                sum(.Object@Tau[,ell]*t(RSS))
              },simplify=FALSE,USE.NAMES=TRUE)
              
              varE.C =  sapply(1:CCO@nbclust,FUN = function(ell){
                RSS = (Wc-.Object@Blup[[ell]]$Nu-.Object@Alpha[[ell]])^2 + .Object@varE/(1+lambda.Nu)
                sum(.Object@Tau[,ell]*RSS)
              },simplify=FALSE,USE.NAMES=TRUE)
              
              .Object@varE = (Reduce("+",varE.D) + Reduce("+",varE.C)) / (n*M)
              
              Gamma2.Nu = sapply(1:CCO@nbclust,FUN = function(ell){
                sum(.Object@Tau[,ell]*.Object@Blup[[ell]]$Nu^2 + .Object@varE /(1+lambda.Nu))
              },simplify=FALSE,USE.NAMES=TRUE)
              .Object@Gamma2.Nu    = Reduce("+",Gamma2.Nu)/n

              Gamma2.Theta = sapply(1:CCO@nbclust,FUN = function(ell){
                Z = (t(.Object@Blup[[ell]]$Theta^2) + .Object@varE /(1+lambda.Theta/Hcoef) ) / Hcoef
                return(t(Z)*.Object@Tau[,ell])
              },simplify=FALSE,USE.NAMES=TRUE)
              
              if (CCO@Gamma2.structure=="scale.location"){
                .Object@Gamma2.Theta = apply(Reduce("+",Gamma2.Theta),2,sum)/n
              }
              if (CCO@Gamma2.structure=="constant"){                
                .Object@Gamma2.Theta = sum(Reduce("+",Gamma2.Theta))/(n*(M-1))
              }
            }
            .Object@eta = getetaFCMM(.Object,CCD,CCO)
            return(.Object) 
          })


setGeneric( name = "getllkFCMM" ,def = function(.Object,CCD,CCO,eta){standardGeneric("getllkFCMM")})
setMethod(
          f          = "getllkFCMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCD,CCO,eta) {
            
            Hcoef = 2^-(CCD@jklevels$scale*eta)
            
            if (CCO@nbclust==1){
              vardens.C           = .Object@Gamma2.Nu+.Object@varE
              vardens.D           = .Object@Gamma2.Theta*Hcoef+.Object@varE    
              logdens   = sapply(CCD@wavecoef, FUN = function(x){
                sum(dnorm( c(x$C,x$D),c(.Object@Alpha,.Object@Beta) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE))
              })
              llk= sum(logdens)
            } else{        
              if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){
                vardens.C   = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Nu[[ell]]+.Object@varE})
                vardens.D   = sapply(1:CCO@nbclust,FUN=function(ell){.Object@Gamma2.Theta[[ell]]*Hcoef+.Object@varE})
              }    
              if ((CCO@Gamma2.structure=="scale.location") | (CCO@Gamma2.structure=="constant")) {
                vardens.C   = .Object@Gamma2.Nu+.Object@varE
                vardens.D   = .Object@Gamma2.Theta*Hcoef+.Object@varE
              }    
              logdens   = sapply(CCD@wavecoef, FUN = function(x){
                apply( dnorm( c(x$C,x$D),sapply(1:CCO@nbclust,FUN = function(ell){c(.Object@Alpha[[ell]],.Object@Beta[[ell]])}) , sd=sqrt(c(vardens.C,vardens.D)),log=TRUE),2,sum)
              })
              
              tau1    = log(.Object@prop)+logdens
              tau.max = apply(tau1,2,max)
              tau1    = exp(t(tau1)-tau.max)
              llk     = sum(log( apply(tau1,1,sum)) + tau.max)
            }  
            return(llk)
          })


setGeneric( name = "getetaFCMM" ,def = function(.Object,CCD,CCO){standardGeneric("getetaFCMM")})
setMethod(
          f          = "getetaFCMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCD,CCO) {
            x1  = 1+1e-4          
            x3  = 6
            r   = (sqrt(5)-1)/2
            x2  = x1+(1-r)*(x3-x1)
            x4  = x1+r*(x3-x1)
            Jx2 = -getllkFCMM(.Object,CCD,CCO,x2)
            Jx4 = -getllkFCMM(.Object,CCD,CCO,x4)
            while (abs(x4-x2)>10^-4){
              if (Jx2<Jx4){
                x3        = x4
                x4        = x2
                Jx4       = Jx2
                x2        = x1+(1-r)*(x3-x1)
                Jx2       = -getllkFCMM(.Object,CCD,CCO,x2)
              } else {
                if (Jx2>Jx4){
                  x1         = x2
                  x2         = x4
                  Jx2        = Jx4
                  x4         = x1+r*(x3-x1)
                  Jx4       = -getllkFCMM(.Object,CCD,CCO,x4)
                } else {
                  x1        = x2
                  x3        = x4
                  x2        = x1+(1-r)*(x3-x1)
                  x4        = x1+r*(x3-x1)        
                  Jx2       = -getllkFCMM(.Object,CCD,CCO,x2)
                  Jx4       = -getllkFCMM(.Object,CCD,CCO,x4)
                }
                
              }
            }
            etamax=max(round((x2+x4)/2,4),1+1e-4)
            etamax
          })

