setMethod(
          f          = "getSNR.mu",
          signature  = "CClustRes",
          definition =
          function(.Object) {

            if (is.list(.Object@Alpha)){
              M      = length(.Object@Beta[[1]]) + 1
              SNR.mu = sapply(1:length(.Object@Alpha),FUN=function(ell){.Object@prop[ell]*(.Object@Alpha[[ell]]^2 + sum(.Object@Beta[[ell]]^2))})/(.Object@varE*M*M)
            } else {
              M      = length(.Object@Beta) + 1
              SNR.mu = (.Object@Alpha^2 + sum(.Object@Beta^2))/(.Object@varE*M*M)
            }                          
            return(SNR.mu)
            })
