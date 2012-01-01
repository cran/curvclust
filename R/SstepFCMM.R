setMethod(
          f          = "SstepFCMM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {
            n              = length(CCD@wavecoef)            
            M              = length(CCD@wavecoef[[1]]$D) + 1
            clustfreq      = apply(.Object@Tau,2,sum)
            .Object@prop   = as.vector(clustfreq) / n 
            Wd             = t(sapply(CCD@wavecoef,  FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc             = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef          = 2^-(CCD@jklevels[,1]*.Object@eta)
            eps            = 1e-10
            
            E     = t(apply(.Object@Tau,1,FUN=function(x) rmultinom(1,1,x)))  
            Nell  = apply(E,2,sum)
            count = 0
            while( (min(Nell)<eps) & (count<100)){
              E     = t(apply(.Object@Tau,1,FUN=function(x) rmultinom(1,1,x)))
              Nell  = apply(E,2,sum)
              count = count + 1
            }

            if (count==100){
              count = 0
              while( (min(Nell)<eps) & (count<10)){
                count = count + 1
                E     = t(apply(matrix(1/CCO@nbclust,ncol=CCO@nbclust,nrow=n),1,FUN=function(x) rmultinom(1,1,x)))   #n*L
                Nell  = apply(E,2,sum)
              }                
            }
            
            .Object@Tau = E            
              if( (CCO@Gamma2.structure=="group.scale.location") | (CCO@Gamma2.structure=="group") ){
                .Object@Blup = sapply(1:CCO@nbclust, FUN=function(ell){
                  Theta = t(apply(.Object@Blup[[ell]]$Theta,1,FUN = function(x){rnorm(M-1,x,sqrt(Hcoef*.Object@varBlup.Theta[[ell]]))}))
                  Nu    = t(t(apply(.Object@Blup[[ell]]$Nu,1,FUN = function(x){rnorm(1,x,sqrt(.Object@varBlup.Nu[[ell]]))})))
                  return(list(Theta=Theta,Nu=Nu))
                },simplify=FALSE,USE.NAMES=TRUE)
              }
              if ((CCO@Gamma2.structure=="scale.location") | (CCO@Gamma2.structure=="constant")) {
                .Object@Blup = sapply(1:CCO@nbclust, FUN=function(ell){
                  Theta = t(apply(.Object@Blup[[ell]]$Theta,1,FUN = function(x){rnorm(M-1,x,sqrt(Hcoef*.Object@varBlup.Theta))}))
                  Nu    = t(t(apply(.Object@Blup[[ell]]$Nu,1,FUN = function(x){rnorm(1,x,sqrt(.Object@varBlup.Nu))})))
                  return(list(Theta=Theta,Nu=Nu))
                },simplify=FALSE,USE.NAMES=TRUE)
            }            
            return(.Object)
          }
          )

