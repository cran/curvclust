setMethod(
          f          = "getwr.mu",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {
            
            M     = length(CCD@wavecoef[[1]]$D) + 1           
            Jmin  = floor(log2(CCD@lengthsignal))
            needs = ceiling(CCD@lengthsignal/2^Jmin)*2^Jmin-CCD@lengthsignal
			M2J   = CCD@lengthsignal+needs
            
            if (CCO@nbclust>1){
              
              mu = lapply(1:CCO@nbclust,FUN = function(ell){
                Wmu.ell = wd(rep(0,M2J),filter.number=CCD@filter.number,family="DaubExPhase")                                
                Wmu.ell = putC.wd(Wmu.ell,level=0,.Object@Alpha[[ell]])
                for (j in unique((CCD@jklevels$scale))){
                  new.coef = rep(0,2^j)
                  pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                  pos4coef = which(CCD@jklevels$scale==j)
                  new.coef[pos4fun] = (.Object@Beta[[ell]][pos4coef])
                  Wmu.ell  = putD.wd(Wmu.ell,level=j,new.coef)
                }
                return(wr(Wmu.ell,filter.number=CCD@filter.number,family="DaubExPhase")[1:CCD@lengthsignal])
              })
            } else {
              
              Wmu = wd(rep(0,M2J),filter.number=CCD@filter.number,family="DaubExPhase")
              Wmu = putC.wd(Wmu,level=0,.Object@Alpha)
              for (j in unique((CCD@jklevels$scale))){
                new.coef = rep(0,2^j)
                pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                pos4coef = which(CCD@jklevels$scale==j)
                new.coef[pos4fun] = (.Object@Beta[pos4coef])
                Wmu  = putD.wd(Wmu,level=j,new.coef)
              }
              mu = wr(Wmu,filter.number=CCD@filter.number,family="DaubExPhase")[1:CCD@lengthsignal]
            }                          
            return(mu)
          })
