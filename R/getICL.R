setMethod(
          f          = "getICL",
          signature  = "CClustRes",
          definition =
          function(.Object,CCD) {
            n             = length(CCD@wavecoef)
            M             = length(CCD@wavecoef[[1]]$D) + 1
            Wd            = t(sapply(CCD@wavecoef,FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc            = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef         = 2^-(CCD@jklevels[,1]*.Object@eta)
            # get nbclust
            if (is.null(.Object@Tau)){nbclust = 1} else {nbclust = dim(.Object@Tau)[2]}
            # get Gamma2.structure
            if (is.list(.Object@Gamma2.Theta)){
              if (length(.Object@Gamma2.Theta[[1]])==1){Gamma2.structure = "group"} else {Gamma2.structure = "group.scale.location"}
            } else {
              if (length(.Object@Gamma2.Theta)==1){Gamma2.structure = "constant"} else {Gamma2.structure = "scale.location"}
            }

            if (nbclust==1){
              Entropy = 0
              if ( (Gamma2.structure=="group.scale.location") | (Gamma2.structure=="scale.location") ) {
                log.RSS.Theta = sum(log( n*.Object@Gamma2.Theta))
                log.RSS.Nu    = log( n*.Object@Gamma2.Nu)
                log.gamma     = -2*M*lgamma(0.5*n)/n                
              } else if ( (Gamma2.structure=="group") | (Gamma2.structure=="constant") ){                                
                log.RSS.Theta = (M-1)*log( n*(M-1)*.Object@Gamma2.Theta)
                log.RSS.Nu    = log( n*.Object@Gamma2.Nu)
                log.gamma     = -2*(lgamma(0.5*n)+lgamma(0.5*n*(M-1)))/n                
              }
            } else {              
              Entropy = -2*sum(.Object@prop*log(.Object@prop))
              Nell    = .Object@prop*n
	      if (min(Nell)>1e-6){
                if (Gamma2.structure=="group.scale.location") {                  
                  log.RSS.Theta = sum(.Object@prop*apply(log(sapply(1:nbclust,FUN=function(ell){.Object@Gamma2.Theta[[ell]]*Nell[ell]})),2,sum))
                  log.RSS.Nu    = sum(.Object@prop*(log(sapply(1:nbclust,FUN=function(ell){.Object@Gamma2.Nu[[ell]]*Nell[ell]}))))
                  log.gamma     = -2*M*sum(lgamma(0.5*Nell))/n                  
              	} else if (Gamma2.structure=="group"){                  
                  log.RSS.Theta = (M-1)*sum(.Object@prop * log(sapply(1:nbclust,FUN = function(ell){.Object@Gamma2.Theta[[ell]]*Nell[ell]*(M-1)})))
                  log.RSS.Nu    = sum(.Object@prop * log(sapply(1:nbclust,FUN = function(ell){.Object@Gamma2.Nu[[ell]]*Nell[ell]})))
                  log.gamma     = -2*sum(lgamma(0.5*Nell) + lgamma(0.5*Nell*(M-1)))/n                  
              	} else if (Gamma2.structure=="scale.location"){                  
                  log.RSS.Theta = sum(log( n*.Object@Gamma2.Theta))
                  log.RSS.Nu    = log( n*.Object@Gamma2.Nu)
                  log.gamma     = -2*M*lgamma(0.5*n)/n                  
                } else if (Gamma2.structure=="constant"){
                  log.RSS.Theta = (M-1)*log( n*(M-1)*.Object@Gamma2.Theta)
                  log.RSS.Nu    = log( n*.Object@Gamma2.Nu)
                  log.gamma     = -2*(lgamma(0.5*n)+lgamma(0.5*n*(M-1)))/n
                }
              } else {
                log.RSS.Theta = NA
                log.RSS.Nu    = NA
                log.gamma     = NA
              }             
            }
            penL    = (M+1)*nbclust*log(n)/n
            ICL  = M*log(n*M*.Object@varE) + log.RSS.Nu + log.RSS.Theta + log.gamma + Entropy + penL                        
            return(ICL)          
          })




