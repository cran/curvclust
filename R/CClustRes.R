setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("listOrNumericOrNULL", c("list","NULL","numeric"))
setClassUnion("matrixOrNULL", c("matrix","NULL"))
setClassUnion("numericOrNULL", c("numeric","NULL"))


setClass(Class = "CClustRes",
         representation(
                        prop          = "numericOrNULL",
                        Beta          = "listOrNumericOrNULL",
                        Alpha         = "listOrNumericOrNULL",
                        Blup          = "listOrNULL",
                        varBlup.Nu    = "listOrNumericOrNULL",
                        varBlup.Theta = "listOrNumericOrNULL",
                        varE          = "numeric",  
                        Gamma2.Theta  = "listOrNumericOrNULL",
                        Gamma2.Nu     = "listOrNumericOrNULL",
                        Tau           = "matrixOrNULL",
                        eta           = "numericOrNULL",
                        loglik        = "numeric"
                        ))



setMethod(
          f          = "initialize",
          signature  = "CClustRes",
          definition =          
          function(.Object,CCO,CCD) {

            n              = length(CCD@wavecoef)            
            M              = length(CCD@wavecoef[[1]]$D) + 1
            .Object@varE   = (mean(CCD@seMAD))^2
            .Object@loglik = 0
            .Object@eta    = 2
            gamma2.init    = (mean(CCD@seMAD))^2
            
            if (CCO@nbclust > 1){
              
              .Object@Tau    = matrix(1/CCO@nbclust,ncol=CCO@nbclust,nrow=n)
              .Object@Beta   = lapply(1:CCO@nbclust,FUN = function(ell){rep(0,M-1)})
              .Object@Alpha  = lapply(1:CCO@nbclust,FUN = function(ell){0})
              .Object@prop   = rep(1/CCO@nbclust,CCO@nbclust)
              
              if (CCO@Gamma2.structure=="group.scale.location"){
                .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                .Object@Gamma2.Theta  = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,M-1)})
                .Object@Gamma2.Nu     = lapply(1:CCO@nbclust,FUN = function(ell){gamma2.init})
                .Object@varBlup.Theta = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,M-1)})
                .Object@varBlup.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){gamma2.init})
                .Object@eta           = 1+1e-4
              }
              
              if (CCO@Gamma2.structure=="group"){
                .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                .Object@Gamma2.Theta  = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                .Object@Gamma2.Nu     = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                .Object@varBlup.Theta = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                .Object@varBlup.Nu    = lapply(1:CCO@nbclust,FUN = function(ell){rep(gamma2.init,1)})
                .Object@eta           = 1+1e-4
              }
              
              if (CCO@Gamma2.structure=="scale.location"){
                .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                .Object@Gamma2.Theta  = rep(gamma2.init,M-1)
                .Object@Gamma2.Nu     = gamma2.init
                .Object@varBlup.Theta = rep(gamma2.init,M-1)
                .Object@varBlup.Nu    = gamma2.init
                .Object@eta           = 1+1e-4
              }            
              
              if (CCO@Gamma2.structure=="constant"){
                .Object@Blup          = lapply(1:CCO@nbclust,FUN = function(ell){list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))})
                .Object@Gamma2.Theta  = gamma2.init
                .Object@Gamma2.Nu     = gamma2.init
                .Object@varBlup.Theta = gamma2.init
                .Object@varBlup.Nu    = gamma2.init
                .Object@eta           = 1+1e-4
              }
              if (CCO@Gamma2.structure=="none"){
                .Object@Blup          = NULL
                .Object@Gamma2.Theta  = NULL
                .Object@Gamma2.Nu     = NULL
                .Object@varBlup.Theta = NULL
                .Object@varBlup.Nu    = NULL
                .Object@eta           = NULL
              }
              
            } else {
              
              .Object@Tau    = NULL
              .Object@Beta   = rep(0,M-1)
              .Object@Alpha  = 0
              .Object@prop   = NULL
              
              if ( (CCO@Gamma2.structure=="scale.location") ){
                .Object@Blup          = list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))
                .Object@Gamma2.Theta  = rep(gamma2.init,M-1)
                .Object@Gamma2.Nu     = gamma2.init
                .Object@varBlup.Theta = rep(gamma2.init,M-1)
                .Object@varBlup.Nu    = gamma2.init
                .Object@eta           = 1+1e-4
              }
              
              if ( (CCO@Gamma2.structure=="constant") ){
                .Object@Blup          = list(Theta = matrix(0,ncol=M-1,nrow=n), Nu=matrix(0,ncol=1,nrow=n))
                .Object@Gamma2.Theta  = gamma2.init
                .Object@Gamma2.Nu     = gamma2.init
                .Object@varBlup.Theta = gamma2.init
                .Object@varBlup.Nu    = gamma2.init
                .Object@eta           = 1+1e-4
              }
              if (CCO@Gamma2.structure=="none"){
                .Object@Blup          = NULL
                .Object@Gamma2.Theta  = NULL
                .Object@Gamma2.Nu     = NULL
                .Object@varBlup.Theta = NULL
                .Object@varBlup.Nu    = NULL
                .Object@eta           = NULL
              }
            }
            return(.Object)
          }
          )


setGeneric( name = "SstepFMM"      ,def = function(.Object,CCO,CCD){standardGeneric("SstepFMM")})
setGeneric( name = "EstepFMM"      ,def = function(.Object,CCO,CCD){standardGeneric("EstepFMM")})
setGeneric( name = "MstepFMM"      ,def = function(.Object,CCO,CCD){standardGeneric("MstepFMM")})
setGeneric( name = "SstepFCMM"     ,def = function(.Object,CCO,CCD){standardGeneric("SstepFCMM")})
setGeneric( name = "EstepFCMM"     ,def = function(.Object,CCO,CCD){standardGeneric("EstepFCMM")})
setGeneric( name = "MstepFCMM"     ,def = function(.Object,CCO,CCD){standardGeneric("MstepFCMM")})
setGeneric( name = "SstepFCM"      ,def = function(.Object,CCO,CCD){standardGeneric("SstepFCM")})
setGeneric( name = "EstepFCM"      ,def = function(.Object,CCO,CCD){standardGeneric("EstepFCM")})
setGeneric( name = "MstepFCM"      ,def = function(.Object,CCO,CCD){standardGeneric("MstepFCM")})
setGeneric( name = "gettheta"      ,def = function(.Object,CCO){standardGeneric("gettheta")})
setGeneric( name = "getBIC"        ,def = function(.Object,CCD){standardGeneric("getBIC")})
setGeneric( name = "getICL"        ,def = function(.Object,CCD){standardGeneric("getICL")})
setGeneric( name = "getwr.mu"      ,def = function(.Object,CCO,CCD){standardGeneric("getwr.mu")})
setGeneric( name = "getwr.Blup"    ,def = function(.Object,CCO,CCD){standardGeneric("getwr.Blup")})
setGeneric( name = "getwr.varU"    ,def = function(.Object,CCO,CCD){standardGeneric("getwr.varU")})
setGeneric( name = "getwr.varBlup" ,def = function(.Object,CCO,CCD){standardGeneric("getwr.varBlup")})
setGeneric( name = "getlambda.U"   ,def = function(.Object,CCO,CCD){standardGeneric("getlambda.U")})
setGeneric( name = "getSNR.mu"     ,def = function(.Object){standardGeneric("getSNR.mu")})


setMethod(f="[","CClustRes",function(x,i,j,drop){
			switch(EXPR=i,
					"prop"         ={return(x@prop)},
					"Beta"         ={return(x@Beta)},
					"Alpha"        ={return(x@Alpha)},
					"Blup"         ={return(x@Blup)},
					"varBlup.Nu"   ={return(x@varBlup.Nu)},
					"varBlup.Theta"={return(x@varBlup.Theta)},
					"varE"         ={return(x@varE)},
					"Gamma2.Theta" ={return(x@Gamma2.Theta)},
					"Gamma2.Nu"    ={return(x@Gamma2.Nu)},
					"Tau"          ={return(x@Tau)},
					"eta"          ={return(x@eta)},
					"loglik"       ={return(x@loglik)},
					stop("This argument does not exist!"))})


setReplaceMethod(f="[","CClustRes",function(x,i,j,value){
			switch(EXPR=i,
					"prop"         ={x@prop<-value},
					"Beta"         ={x@Beta<-value},
					"Alpha"        ={x@Alpha<-value},
					"Blup"         ={x@Blup<-value},
					"varBlup.Nu"   ={x@varBlup.Nu<-value},
					"varBlup.Theta"={x@varBlup.Theta<-value},
					"varE"         ={x@varE<-value},
					"Gamma2.Theta" ={x@Gamma2.Theta<-value},
					"Gamma2.Nu"    ={x@Gamma2.Nu<-value},
					"Tau"          ={x@Tau<-value},
					"eta"          ={x@eta<-value},
					"loglik"       ={x@loglik<-value},
					stop("This argument does not exist!"))
			validObject(x)
			return(x)}) 




setMethod(
		f = "summary",
		signature = "CClustRes",
		definition = function(object){			
			if (is.null(object@Tau)){
				nbclust = 1
			} else {
				nbclust = dim(object@Tau)[2]
			}
			cat("Number of clusters ",nbclust,"\n")
			if (is.null(object@Tau)){
				csize = 1
			} else {
				csize = object@prop
			cat("cluster size ",csize,"\n")
			cat("labels ","\n")
			print(apply(object@Tau,1,which.max))
			}

		}          
)
