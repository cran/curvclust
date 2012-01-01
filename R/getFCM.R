setMethod(
          f          = "getFCM",
          signature  = "CClustData",
          definition =
          function(.Object,CCO) {

            if (CCO@nbclust==1){
              stop("if CCO@nbclust==1 use getFMM instead")
            }
            if (CCO@Gamma2.structure!="none"){
              stop("if CCO@Gamma2.structure!=\"none\" use getFCMM instead")
            }
            
            cat("[ initialization with",CCO@burn,"burns using",CCO@init,"]\n") 
            llik      = c()
            loglikbm1 = -Inf
            CCR       = new("CClustRes",CCO,.Object)                        
            CCRu      = new("CClustRes",CCO,.Object)
            CCRum1    = new("CClustRes",CCO,.Object)

            if (CCO@init=="SEM"){
              b = 0
              while ( (b < CCO@burn)){
                b   = b+1
                CCR = SstepFCM(CCR,CCO,.Object)
                CCR = MstepFCM(CCR,CCO,.Object)
                CCR = EstepFCM(CCR,CCO,.Object)
                if ( (CCR@loglik> loglikbm1) & (CCR@varE>1e-10) & (min(CCR@prop)>1e-6)) {
                  CCRum1    = CCRu
                  CCRu      = CCR
                  loglikbm1 = CCR@loglik
                }
                llik = c(llik,CCR@loglik)              
              }
              if (min(CCR@prop)<1e-10){stop("[SEM initializations resulted in empty clusters]\n")}
              if (CCR@varE<1e-10){stop("[SEM initializations resulted in null variance]\n")}
              
            } else if (CCO@init=="rEM"){
              dv = 0
              for (nbseeds in c(1:CCO@burn)){
                CCR  = new("CClustRes",CCO,.Object)
                CCR  = SstepFCM(CCR,CCO,.Object)
                i    = 0
                while ( (i < 10) & (min(CCR@prop)>1e-10) & (CCR@varE>1e-10)){
                  i   = i+1
                  CCR = MstepFCM(CCR,CCO,.Object)
                  CCR = EstepFCM(CCR,CCO,.Object)
                }
                if (min(CCR@prop)<1e-10){CCR@loglik=-Inf; dv=dv+1}
                if (CCR@loglik> loglikbm1){
                  CCRum1    = CCRu
                  CCRu      = CCR
                  loglikbm1 = CCR@loglik
                }
                llik = c(llik,CCR@loglik)
              }
              if (dv==CCO@burn){stop("[All rEM initializations resulted in empty clusters]\n")}
            }

            cat("[EM algorithm]\n")            
            CCR           = CCRu            
            theta.h       = gettheta(CCRu,CCO)
            theta.hm1     = gettheta(CCRum1,CCO)
            theta.dot.hm2 = gettheta(CCRu,CCO)
            epsilon       = Inf
            
            while ( (epsilon>CCO@eps) & (min(CCR@prop)>1e-10)  & (CCR@varE>1e-10)){
              
              CCR           = MstepFCM(CCR,CCO,.Object)
              CCR           = EstepFCM(CCR,CCO,.Object)
              
              theta.hp1     = gettheta(CCR,CCO)              
              theta.dot.hm1 = theta.h + invnorm( invnorm(theta.hm1-theta.h) + invnorm( theta.hp1- theta.h) )
              epsilon       = sum( (theta.dot.hm1-theta.dot.hm2)^2 / (theta.dot.hm1)^2 )
              theta.hm1     = theta.h
              theta.h       = theta.hp1
              theta.dot.hm2 = theta.dot.hm1
              llik          = c(llik,CCR@loglik)
            }

            if (min(CCR@prop)<1e-10){cat("[EM stopped because of empty clusters]")}            
            if (CCO@loglikplot){plot(llik,type="l");abline(v=CCO@burn)}
            if (CCR@varE <1e-10){cat("[EM stopped because of null variance]")}
            
            return(CCR)
          })

invnorm <- function(x){
  y = x/sum(x^2)
  y[is.na(y)] = 0
  return(y)
}

