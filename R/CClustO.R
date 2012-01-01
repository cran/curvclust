setClassUnion("charOrNULL", c("character","NULL"))


setClass(Class = "CClustO",representation(
				nbclust            = "numeric",
				Gamma2.structure   = "charOrNULL",
				burn               = "numeric",
				eps                = "numeric",
				init               = "character",
				loglikplot         = "logical"
		),
		prototype(nbclust = 2,
				Gamma2.structure = "constant",
				burn             = 100,
				eps              = 0.1,
				init             = "SEM",
				loglikplot       = FALSE
		),         
		validity=function(object){
		stopifnot(object@burn>=0)
		stopifnot(object@eps>=0)
		stopifnot(object@Gamma2.structure %in% c("none","constant", "group", "scale.location", "group.scale.location"))
		stopifnot(object@init %in% c("rEM","SEM"))
			return(TRUE)         
		})


setMethod(f="[","CClustO",function(x,i,j,drop){
			switch(EXPR=i,
					"nbclust"         ={return(x@nbclust)},
					"Gamma2.structure"={return(x@Gamma2.structure)},
					"burn"            ={return(x@burn)},
					"eps"             ={return(x@eps)},
					"loglikplot"      ={return(x@loglikplot)},
					"init"            ={return(x@init)},
					stop("This argument does not exist!"))})

setReplaceMethod(f="[","CClustO",function(x,i,j,value){
			switch(EXPR=i,
					"nbclust"         ={x@nbclust<-value},
					"Gamma2.structure"={x@Gamma2.structure<-value},
					"burn"            ={x@burn<-value},
					"eps"             ={x@eps<-value},
					"loglikplot"      ={x@loglikplot<-value},
					"init"            ={x@init<-value},
					stop("This argument does not exist!"))
			validObject(x)
			return(x)})



