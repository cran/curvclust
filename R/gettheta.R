setMethod(
          f          = "gettheta",
          signature  = "CClustRes",
          definition =
          function(.Object,CCO) {
            
            if ((CCO@nbclust>1) & (CCO@Gamma2.structure!="none")){              
              theta = c(.Object@prop,unlist(.Object@Beta),unlist(.Object@Alpha),.Object@varE,unlist(.Object@Gamma2.Theta),unlist(.Object@Gamma2.Nu),.Object@eta)              
            }
            if ((CCO@nbclust==1) & (CCO@Gamma2.structure!="none")){
              theta = c(unlist(.Object@Beta),unlist(.Object@Alpha),.Object@varE,unlist(.Object@Gamma2.Theta),unlist(.Object@Gamma2.Nu),.Object@eta)              
            }
            if ((CCO@nbclust>1) & (CCO@Gamma2.structure=="none")){
              theta = c(.Object@prop,unlist(.Object@Beta),unlist(.Object@Alpha),.Object@varE,.Object@eta)               
            }
            return(theta)
          }
          )

