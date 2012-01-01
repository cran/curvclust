setMethod(
          f          = "getwr.Blup",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {

			M     = length(CCD@wavecoef[[1]]$D) + 1           
            Jmin  = floor(log2(CCD@lengthsignal))
            needs = ceiling(CCD@lengthsignal/2^Jmin)*2^Jmin-CCD@lengthsignal
			M2J   = CCD@lengthsignal+needs
			
			            
            if (CCO@nbclust>1){
                          
              Blup = lapply(.Object@Blup,FUN = function(B){
                U.ell = sapply(1:length(names(CCD@wavecoef)),FUN =function(i){
                  WUi.ell = wd(rep(0,M2J),filter.number=CCD@filter.number,family="DaubExPhase")
                  WUi.ell = putC.wd(WUi.ell,level=0,B$Nu[i,1])
                  for (j in unique((CCD@jklevels$scale))){
                    new.coef = rep(0,2^j)
                    pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                    pos4coef = which(CCD@jklevels$scale==j)
                    new.coef[pos4fun] = (B$Theta[i,pos4coef])
                    WUi.ell  = putD.wd(WUi.ell,level=j,new.coef)
                  }
                  return(wr(WUi.ell,filter.number=CCD@filter.number,family="DaubExPhase"))
                })
                colnames(U.ell)=names(CCD@wavecoef)
                return(t(U.ell))
              })

            } else {             
              
              Blup = lapply(1:length(names(CCD@wavecoef)),FUN =function(i){
                WUi = wd(rep(0,M2J),filter.number=CCD@filter.number,family="DaubExPhase")
                WUi = putC.wd(WUi,level=0,.Object@Blup$Nu[i,1])
                for (j in unique((CCD@jklevels$scale))){
                  new.coef = rep(0,2^j)
                  pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                  pos4coef = which(CCD@jklevels$scale==j)
                  new.coef[pos4fun] = (.Object@Blup$Theta[i,pos4coef])
                  WUi  = putD.wd(WUi,level=j,new.coef)
                }
                return(wr(WUi,filter.number=CCD@filter.number,family="DaubExPhase"))
              })
              names(Blup) = names(CCD@wavecoef)
            }                          
            return(Blup)
            })
