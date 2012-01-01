setMethod(
          f          = "getlambda.U",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {
            
            if (is.null(.Object@Blup)){
            	lambda.U = NULL
            } else if (is.list(.Object@Alpha)){
              M      = length(.Object@Beta[[1]]) + 1
              if (CCO@Gamma2.structure=="constant"){                
                lambda.U = M*.Object@varE/(.Object@Gamma2.Nu + .Object@Gamma2.Theta/(1-2^(1-.Object@eta)))                
              } else if (CCO@Gamma2.structure=="group"){
                lambda.U = lapply(1:length(.Object@Alpha),FUN=function(ell){M*.Object@varE/(.Object@Gamma2.Nu[[ell]] + .Object@Gamma2.Theta[[ell]]/(1-2^(1-.Object@eta)))})
              } else if (CCO@Gamma2.structure=="scale.location"){
                lambda.U = M*.Object@varE/(.Object@Gamma2.Nu+sum(2^(-.Object@eta*CCD@jklevels$scale)*.Object@Gamma2.Theta))
              } else if (CCO@Gamma2.structure=="group.scale.location"){
                lambda.U = lapply(1:length(.Object@Alpha),FUN=function(ell){M*.Object@varE/(.Object@Gamma2.Nu[[ell]] + sum(2^(-.Object@eta*CCD@jklevels$scale)*.Object@Gamma2.Theta[[ell]]))})
              }
            } else {
              M      = length(.Object@Beta) + 1
              if (CCO@Gamma2.structure=="constant"){
                lambda.U = M*.Object@varE/(.Object@Gamma2.Nu + .Object@Gamma2.Theta/(1-2^(1-.Object@eta)))
              } else if (CCO@Gamma2.structure=="scale.location"){
                lambda.U = M*.Object@varE/(.Object@Gamma2.Nu+sum(2^(-.Object@eta*CCD@jklevels$scale)*.Object@Gamma2.Theta))
              }
            }
            return(lambda.U)
          })          
          
