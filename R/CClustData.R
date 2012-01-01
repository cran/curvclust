setClassUnion("numOrNULL", c("numeric","NULL"))
setClassUnion("matOrNULL", c("matrix","NULL"))

anyMissing = function (x){any(is.na(unlist(x)))}

fillNA     = function(v){v[which(is.na(v))]=median(v[which(!(is.na(v)))]);v}

setClass(Class = "CClustData", representation  (
           wavecoef      = "list",
           jklevels      = "data.frame",
           seMAD         = "numeric",
           lengthsignal  = "numeric",
           filter.number = "numeric"
           )
         )
         
setMethod(
          f          = "initialize",
          signature  = "CClustData",
          definition =          
          function(.Object,Y,filter.number) {

            if(missing(filter.number)){filter.number=8}

            n  = length(Y)
            if(is.null(names(Y))){
              names(Y) = paste("ID",1:n,sep="")
            }

            mi = sapply(Y,FUN = function(y){length(y)})
            if (length(unique(mi))>1){
              stop("All curves should have the same number of recorded positions")
            }
            if(any(is.na(unlist(Y)))){stop("please remove missing values")}

            Minit = length(Y[[1]])
            .Object@lengthsignal = length(Y[[1]])
            
            Jmin  = floor(log2(Minit))
            needs = ceiling(Minit/2^Jmin)*2^Jmin-Minit

            if (needs>0){
              Y = sapply(Y,FUN = function(y){y = c(y,rep(0,needs)); y},simplify=FALSE,USE.NAMES=TRUE)
            }
          
            M      = length(Y[[1]])
            J      = log2(M)-1
            
            wavdec = sapply(Y,FUN=function(y){wd(y,filter.number=filter.number,family="DaubExPhase")},simplify=F,USE.NAMES=TRUE)
            D    = sapply(wavdec, FUN=function(x){x$D},simplify=F,USE.NAMES=TRUE)
            C    = unlist(sapply(wavdec, FUN=function(x){accessC(x,level=0)},simplify=F,USE.NAMES=TRUE))
            
            jcoef              = rep(J:0,2^(J:0))
            kcoef              = unlist(lapply(J:0,FUN=function(x) 0:(2^x-1)))
            jklevels           = data.frame(scale=jcoef,location=kcoef)

            seMAD = sapply(Y, FUN = function(y){
              u         = floor(log2(length(y)))
              ywd       = wd(y,filter.number=filter.number,family="DaubExPhase")
              finecoef  = accessD(ywd,lev=(u-1))
              sigma     = mad(finecoef) 
            })

            .Object@seMAD = seMAD
       
            if(needs>0){
              pos2keep  = apply(Reduce("cbind",D),1,FUN=function(x){sum(x==0)})!=n
              jklevels  = jklevels[pos2keep,]
            } else {
              pos2keep = rep(TRUE,M-1)
            }

            .Object@wavecoef = sapply(names(D),FUN = function(ID){
              list(C=C[[ID]],D=D[[ID]][pos2keep])
            },simplify=F,USE.NAMES=TRUE)

            .Object@jklevels      = jklevels            
            .Object@filter.number = filter.number
            
            return(.Object)
          }
          )



setGeneric( name = "getFMM"   ,def = function(.Object,CCO){standardGeneric("getFMM")})
setGeneric( name = "getFCMM"  ,def = function(.Object,CCO){standardGeneric("getFCMM")})
setGeneric( name = "getFCM"   ,def = function(.Object,CCO){standardGeneric("getFCM")})
setGeneric( name = "getCM"    ,def = function(.Object,CCO){standardGeneric("getCM")})
setGeneric( name = "getUnionCoef"   ,def = function(.Object){standardGeneric("getUnionCoef")})
setGeneric( name = "getwr.Y"   ,def = function(.Object){standardGeneric("getwr.Y")})
setGeneric( name = "getwr"   ,def = function(.Object){standardGeneric("getwr")})

setMethod(f="[","CClustData",function(x,i,j,drop){
			switch(EXPR=i,
					"wavecoef"       ={return(x@wavecoef)},
					"jklevels"      ={return(x@jklevels)},
					"seMAD"         ={return(x@seMAD)},
					"lengthsignal"  ={return(x@lengthsignal)}, 
					"filter.number" ={return(x@filter.number)},
					stop("This argument does not exist!"))})         

setReplaceMethod(f="[","CClustData",function(x,i,j,value){
			switch(EXPR=i,
					"wavecoef"       ={x@wavecoef<-value},
					"jklevels"      ={x@jklevels<-value},
					"seMAD"         ={x@seMAD<-value},
					"lengthsignal"  ={x@lengthsignal<-value}, 
					"filter.number" ={x@filter.number<-value},
					stop("This argument does not exist!"))
			validObject(x)
			return(x)})



setMethod(
		f = "summary",
		signature = "CClustData",
		definition = function(object){
		cat("Number of observations ", length(object@wavecoef),"\n")		
		cat("Number of positions ", length(object@wavecoef[[1]]$D) + 1,"\n")		
		cat("Filter Number ",object@filter.number,"\n")
		cat("positions\n")
		for (j in unique(object@jklevels$scale)){
			cat("scale j=",j,"\n")
			cat("positions k(j)=",object@jklevels[object@jklevels$scale==j,]$location,"\n")			
		}		
		}          
)

setMethod(
		f = "show",
		signature = "CClustData",
		definition = function(object){
			cat("***************************\n")			
			cat("Data are in the list format\n")
			cat("scaling coefficients\n")
			cat("CCD[\"wavecoef\"][[1]]$C\n")
			print(object@wavecoef[[1]]$C)
			cat("wavelet coefficients\n")
			cat("CCD[\"wavecoef\"][[1]]$D[1:10]\n")
			print(object@wavecoef[[1]]$D[1:10])
		}          
)
