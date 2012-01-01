setMethod(
          f          = "getwr.Y",
          signature  = "CClustData",
          definition =
          function(.Object) {
            
            M     = length(.Object@wavecoef[[1]]$D) + 1           
            Jmin  = floor(log2(.Object@lengthsignal))
            needs = ceiling(.Object@lengthsignal/2^Jmin)*2^Jmin-.Object@lengthsignal
			M2J   = .Object@lengthsignal+needs
                                                        
            Y = lapply(.Object@wavecoef,FUN = function(x){
              WY = wd(rep(0,M2J),filter.number=.Object@filter.number,family="DaubExPhase")
              WY = putC.wd(WY,level=0,x$C)
              for (j in unique((.Object@jklevels$scale))){
                new.coef = rep(0,2^j)
                pos4fun  = .Object@jklevels[.Object@jklevels$scale==j,2]+1
                pos4coef = which(.Object@jklevels$scale==j)
                new.coef[pos4fun] = (x$D[pos4coef])
                WY = putD.wd(WY,level=j,new.coef)
              }
              return(wr(WY,filter.number=.Object@filter.number,family="DaubExPhase")[1:.Object@lengthsignal])
            })
            return(Y)
          })
            
