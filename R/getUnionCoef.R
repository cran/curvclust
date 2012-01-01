setMethod(
          f          = "getUnionCoef",
          signature  = "CClustData",
          definition =
          function(.Object) {

            M  = dim(.Object@jklevels)[1]+1
            n  = length(names(.Object@wavecoef))
            se = mean(.Object@seMAD)            
            zeros = sapply(.Object@wavecoef,FUN = function(x){
              coef = pmax(unlist(x)-se*sqrt(2*log(M-1)),0)
              return(coef)
            })
            union    = (apply(zeros,1,sum)!=0)
            union[1] = TRUE
            newcoef = lapply(.Object@wavecoef,FUN = function(x){
              z = c(x$C,x$D)[union]
              C = z[1]
              D = z[-1]
              list(C=C,D=D)
              list(C=C,D=D)
            })            
            .Object@wavecoef = newcoef
            .Object@jklevels = .Object@jklevels[union[-1],]
            return(.Object)
          }          
          )
