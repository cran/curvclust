setMethod(
          f          = "MstepFCM",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {

            tmp           = order(unlist(.Object@Alpha))
            .Object@Beta  = .Object@Beta[tmp]
            .Object@Alpha = .Object@Alpha[tmp]
            .Object@prop  = .Object@prop[tmp]
            .Object@Tau   = .Object@Tau[,tmp]
            
            n             = length(CCD@wavecoef)
            M             = length(CCD@wavecoef[[1]]$D) + 1
            Wd            = t(sapply(CCD@wavecoef,FUN=function(x){x$D},USE.NAMES=TRUE))
            Wc            = t(t(sapply(CCD@wavecoef,FUN=function(x){x$C},USE.NAMES=TRUE)))
            Hcoef         = 2^-(CCD@jklevels[,1]*.Object@eta)

            clustfreq     = apply(.Object@Tau,2,sum)
            .Object@prop  = as.vector(clustfreq) / n
            
            .Object@Beta  = lapply(1:CCO@nbclust,FUN = function(ell){
              apply(.Object@Tau[,ell]*Wd,2,sum)/clustfreq[ell]
            })
            .Object@Alpha = lapply(1:CCO@nbclust,FUN = function(ell){
              sum(.Object@Tau[,ell]*Wc)/clustfreq[ell]
            })

            varE.D =  sapply(1:CCO@nbclust,FUN = function(ell){
              apply( .Object@Tau[,ell]* t(apply(Wd,1,FUN = function(x){x-.Object@Beta[[ell]]})^2) ,2,sum)
            },simplify=FALSE,USE.NAMES=TRUE)
            
            varE.C =  sapply(1:CCO@nbclust,FUN = function(ell){
              sum(.Object@Tau[,ell]*( (Wc - .Object@Alpha[[ell]])^2))
            },simplify=FALSE,USE.NAMES=TRUE)
            
            .Object@varE = sum( (Reduce("+",varE.D) + Reduce("+",varE.C))) / (n*M)
            return(.Object) 
          })
