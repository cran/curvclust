setMethod(
          f          = "getBIC",
          signature  = "CClustRes",
          definition =
          function(.Object,CCD) {

            n = length(CCD@wavecoef)            
            M = length(CCD@wavecoef[[1]]$D) + 1            

            if (is.null(.Object@Tau)){
              nbclust = 1
            } else {
              nbclust = dim(.Object@Tau)[2]
            }
            
            if (is.null(.Object@Gamma2.Theta)){
              dimL = (M+1)*nbclust
            } else {
              if (is.list(.Object@Gamma2.Theta)){
                if (length(.Object@Gamma2.Theta[[1]])==1){
                  Gamma2.structure = "group"
                } else {
                  Gamma2.structure = "group.scale.location"                
                }
              } else {
                if (length(.Object@Gamma2.Theta)==1){
                  Gamma2.structure = "constant"
                } else {
                  Gamma2.structure = "scale.location"                
                }
              }
              
              if( (Gamma2.structure=="group.scale.location") ){
                dimL = ( (M+1)*nbclust + M*nbclust )                 
              } else if (Gamma2.structure=="group"){
                dimL = ( (M+1)*nbclust + 2*nbclust )                        
              } else if (Gamma2.structure=="scale.location"){
                dimL = ( (M+1)*nbclust + M )                        
              } else if (Gamma2.structure=="constant"){
                dimL = ( (M+1)*nbclust + 2 )                        
              }
            }            
            
            BIC = .Object@loglik - 0.5*dimL*log(n)
            
            return(BIC) 
          })
