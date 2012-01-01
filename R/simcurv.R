simcurv <- function(mu,prop,SNR,lambda,n,eta,filter.number,Gamma2.structure="constant"){
	
	
	prop    = sort(prop,decreasing=TRUE)
	M       = dim(mu)[1]
	L       = dim(mu)[2]
	J       = log(M)/log(2)
	s2e     = sum(prop*t(t(t(mu)-apply(mu,2,mean))^2))/((M*SNR)^2)*M
	label   = apply(rmultinom(n,1,prop),2,which.max)
	
	i = 0
	while (length(unique(label))!= length(prop) & (i<1000)){
		label = apply(rmultinom(n,1,prop),2,which.max)
		i = i+1
	}
	
	if (length(unique(label))!= length(prop) & (i==1000)){
		stop("[All samplings resulted in empty clusters]\n")
	}
	
	Tau     = matrix(0,ncol=length(prop),nrow=n)
	for (i in c(1:n)){Tau[i,label[i]]=1}
	transfo = list()
	
	for (ell in c(1:L)){
		eval(parse(text=paste("transfo[[",ell,"]]=wd(mu[,",ell,"],filter.number=",filter.number,",family=\"DaubExPhase\")",sep="")))
	}
	
	mu.noise = lapply(label,FUN = function(ell){
				E       = rnorm(M,0,sqrt(s2e))
				epsilon = wd(E,filter.number=filter.number,family="DaubExPhase")
				beta    = wd(rep(0,M),filter.number=filter.number,family="DaubExPhase")
				beta    = putC.wd(beta,level=0,accessC(transfo[[ell]],level=0)+accessC(epsilon,level=0))
				for (jlevel in c(1:J)){    
					beta = putD.wd(beta,level=(jlevel-1),accessD(transfo[[ell]],level=(jlevel-1))+accessD(epsilon,level=(jlevel-1)))
				}
				return(wr(beta,filter.number=filter.number,family="DaubExPhase"))
			})	
	  
	CCO    = new("CClustO")
	CCO@Gamma2.structure = Gamma2.structure

	
	if (Gamma2.structure=="none"){
		CCD       = new("CClustData",lapply(1:n,FUN = function(i){mu.noise[[i]]}),filter.number)
		CCR       = new("CClustRes",CCO,CCD)
		CCR@Tau   = Tau
		CCR@prop  = apply(Tau,2,sum)/n
		
		wavecoef  = apply(mu,2,FUN=function(y){wd(y,filter.number=filter.number,family="DaubExPhase")})
		CCR@Beta  = sapply(wavecoef, FUN=function(x){x$D},simplify=F,USE.NAMES=TRUE)
		CCR@Alpha = sapply(wavecoef, FUN=function(x){accessC(x,level=0)},simplify=F,USE.NAMES=TRUE)
		CCR@varE  = s2e
		  
	} else {
		if (Gamma2.structure=="constant"){
			gamma2  = M*s2e/(lambda*(1+1/(1-2^(1-eta))))
			U = lapply(label,FUN = function(ell){    
						theta = wd(rep(0,M),filter.number=filter.number,family="DaubExPhase")
						theta = putC.wd(theta,level=0,rnorm(1,0,sqrt(median(s2e))))
						for (jlevel in c(1:J)){
							theta = putD.wd(theta,level=(jlevel-1),rnorm(2^(jlevel-1),0,sqrt(2^(-(jlevel-1)*eta)*gamma2)))
						}
						return(wr(theta,filter.number=filter.number,family="DaubExPhase"))
					})
		} else if (Gamma2.structure=="scale.location"){
			gamma2 = rgamma(M, shape = M*s2e/(lambda*(1+1/(1-2^(1-eta))))/2, scale = 2)
			jcoef  = rep((J-1):0,2^((J-1):0))
			U       = lapply(label,FUN = function(ell){    
						theta = wd(rep(0,M),filter.number=filter.number,family="DaubExPhase")
						theta = putC.wd(theta,level=0,rnorm(1,0,sqrt(median(s2e))))
						for (jlevel in c(1:J)){
							theta = putD.wd(theta,level=(jlevel-1),rnorm(2^(jlevel-1),0,sqrt(2^(-(jlevel-1)*eta)*gamma2[which(jcoef==(jlevel-1))])))
						}
						return(wr(theta,filter.number=filter.number,family="DaubExPhase"))
					})
		}
		
		CCD      = new("CClustData",lapply(1:n,FUN = function(i){mu.noise[[i]]+U[[i]]}),filter.number)
		CCR       = new("CClustRes",CCO,CCD)
		CCR@Tau   = Tau
		CCR@prop  = apply(Tau,2,sum)/n
		
		wavecoef  = apply(mu,2,FUN=function(y){wd(y,filter.number=filter.number,family="DaubExPhase")})
		CCR@Beta  = sapply(wavecoef, FUN=function(x){x$D},simplify=F,USE.NAMES=TRUE)
		CCR@Alpha = sapply(wavecoef, FUN=function(x){accessC(x,level=0)},simplify=F,USE.NAMES=TRUE)
		CCR@varE  = s2e
		
		Blup     = lapply(U,FUN=function(u){wd(u,filter.number=filter.number,family="DaubExPhase")})
		ThetaM   = t(sapply(Blup, FUN=function(x){x$D}))
		NuM      = t(t(sapply(Blup, FUN=function(x){accessC(x,level=0)})))
		Theta    = lapply(1:L,FUN = function(ell){a = as.numeric(label==ell);return(ThetaM*a)})
		Nu       = lapply(1:L,FUN = function(ell){a = as.numeric(label==ell);return(NuM*a)}) 
		CCR@Blup = lapply(1:L,FUN=function(ell) list(Theta = Theta[[ell]] ,Nu = Nu[[ell]])) 		
		CCR@Gamma2.Theta = gamma2
		CCR@Gamma2.Nu    = sqrt(median(s2e))			
	}
	
	return(list(CCD=CCD,CCR=CCR))
	
}


