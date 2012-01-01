setMethod(
          f          = "getwr.varU",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO,CCD) {

			M     = length(CCD@wavecoef[[1]]$D) + 1           
            Jmin  = floor(log2(CCD@lengthsignal))
            needs = ceiling(CCD@lengthsignal/2^Jmin)*2^Jmin-CCD@lengthsignal
			M2J   = CCD@lengthsignal+needs
            W     = t(GenW(M2J,filter.number=CCD@filter.number,family="DaubExPhase"))
            
            if (CCO@nbclust>1){
              
              if (CCO@Gamma2.structure=="constant"){
                pos = .Object@Gamma2.Nu
                for (j in c((log2(M2J)-1):0)){
                  new.coef = rep(0,2^j)
                  pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                  pos4coef = which(CCD@jklevels$scale==j)
                  new.coef[pos4fun] = .Object@Gamma2.Theta*2^(-j*.Object@eta)
                  pos      = c(pos,new.coef)
                }                
                varU = (t(W)%*%diag(pos)%*%(W))[1:CCD@lengthsignal,1:CCD@lengthsignal]
                
                
              } else if (CCO@Gamma2.structure=="group"){
                varU = lapply(1:CCO@nbclust,FUN=function(ell){
                  pos  = .Object@Gamma2.Nu[[ell]]                                
                  for (j in c((log2(M2J)-1):0)){
                    new.coef = rep(0,2^j)                    
                    pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                    pos4coef = which(CCD@jklevels$scale==j)
                    new.coef[pos4fun] = .Object@Gamma2.Theta[[ell]]*2^(-j*.Object@eta)
                    pos      = c(pos,new.coef)
                  }
                  return((t(W)%*%diag(pos)%*%(W))[1:CCD@lengthsignal,1:CCD@lengthsignal])
                })
                
              } else if (CCO@Gamma2.structure=="scale.location"){
                pos  = .Object@Gamma2.Nu
                for (j in c((log2(M2J)-1):0)){              
                  new.coef = rep(0,2^j)
                  pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                  pos4coef = which(CCD@jklevels$scale==j)
                  new.coef[pos4fun] = .Object@Gamma2.Theta[pos4coef]*2^(-j*.Object@eta)
                  pos = c(pos,new.coef)
                }        
                varU = (t(W)%*%diag(pos)%*%(W))[1:CCD@lengthsignal,1:CCD@lengthsignal]
                
              } else if (CCO@Gamma2.structure=="group.scale.location"){
                varU = lapply(1:CCO@nbclust,FUN=function(ell){
                  pos  = .Object@Gamma2.Nu[[ell]]
                  for (j in c((log2(M2J)-1):0)){
                    new.coef = rep(0,2^j)                    
                    pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                    pos4coef = which(CCD@jklevels$scale==j)
                    new.coef[pos4fun] = .Object@Gamma2.Theta[[ell]][pos4coef]*2^(-j*.Object@eta)
                    pos      = c(pos,new.coef)
                  }
                  return((t(W)%*%diag(pos)%*%(W))[1:CCD@lengthsignal,1:CCD@lengthsignal])
                })                
              }
            } else {
              if (CCO@Gamma2.structure=="constant"){
                pos = .Object@Gamma2.Nu
                for (j in c((log2(M2J)-1):0)){
                  new.coef = rep(0,2^j)
                  pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                  pos4coef = which(CCD@jklevels$scale==j)
                  new.coef[pos4fun] = .Object@Gamma2.Theta*2^(-j*.Object@eta)
                  pos      = c(pos,new.coef)
                }                
                varU = (t(W)%*%diag(pos)%*%(W))[1:CCD@lengthsignal,1:CCD@lengthsignal]
              } else if (CCO@Gamma2.structure=="scale.location"){
                pos  = .Object@Gamma2.Nu
                for (j in c((log2(M2J)-1):0)){
                  new.coef = rep(0,2^j)
                  pos4fun  = CCD@jklevels[CCD@jklevels$scale==j,2]+1
                  pos4coef = which(CCD@jklevels$scale==j)
                  new.coef[pos4fun] = .Object@Gamma2.Theta[pos4coef]*2^(-j*.Object@eta)
                  pos = c(pos,new.coef)
                }
                varU = (t(W)%*%diag(pos)%*%(W))[1:CCD@lengthsignal,1:CCD@lengthsignal]             
              }
            }
            return(varU)
          })
