isTest <-  function(isObj,kStar=NULL,alpha=0.05){
	#browser()
	w <- weights(isObj)
	N <- length(w)
	if(is.null(kStar)) kStar <- 3*N/5
	u <- sort(w)[N-kStar]
	z <- w-u
	z<-z[z>0]
	n<- length(z)
	
	fLogLikUnrestr <- function(tau){
		csi <-mean(log(1+tau*z))
		(log(csi/tau)+(1+1/csi)*mean(log(1+z*tau)))
		}
	fLogLikUnrestr <-Vectorize(fLogLikUnrestr)
	
	opt <- nlm(f= fLogLikUnrestr,p=1/mean(z))
	logLikUnrestr <- n*opt$minimum
	tau <- opt$estimate
	csi <- mean(log(1+tau*z))
	beta <- csi/tau
	
	print(csi)
	
	if(csi<=0.5) warning("succede!!!")

	t <- sqrt(n/(3*beta^2))*(csi-0.5)
	pW <- pnorm(t,lower.tail=FALSE)
	
	fLogLikRestr <- function(beta){
		 if(beta<=sqrt(.Machine$double.eps))return(Inf)
		 (log(beta)+3*mean(log(1+z/(2*beta))))
		 }
	fLogLikRestr <-Vectorize(fLogLikRestr)
	# plot(fLogLikRestr)
	
	opt <- nlm(f=fLogLikRestr,p=mean(z)/2)	
	logLikRestr <- n*opt$minimum
	betaH0 <- opt$estimate
	
	sCsi <- (4*sum(log(1+z/(2*betaH0)))-6*sum(z/(2*betaH0+z)))/sqrt(2*n)
	pS <- pnorm(sCsi,lower.tail=FALSE)
	
	LR <- 2*(logLikRestr-logLikUnrestr)
	pLRT <- 1-(pchisq(LR,0)+pchisq(LR,1))/2
	
	return(new("isTest",waldT=c("p.value"=pW,"testStat"=t,"logLik"=logLikUnrestr),
	                 scoreT=c("p.value"=pS,"testStat"=sCsi,"logLik"= logLikRestr),
	                 LRT = c("p.value"=pLRT,"testStat"=LR),
	                 para=c("csi"=csi,"beta"=beta,"betaH0"=betaH0),
	                 logLikRestr=fLogLikRestr,
	                 logLikUnRestr= fLogLikUnrestr,
	                 kStar=kStar,
	                 alpha=alpha,
	                 n = n))
	
	}
	
	
	
setClass("isTest",representation=representation(waldT="numeric",
                                                scoreT="numeric",
                                                LRT = "numeric",
                                                para="numeric",
                                                logLikUnRestr="function",
                                                logLikRestr="function",
                                                kStar="numeric",
                                                alpha="numeric",
                                                n = "integer"))


### Minimal description
setMethod("show",signature("isTest"),function(object){
	# cat("Wald t-test:\n")
	print(rbind("t"=object@waldT,"score"=object@scoreT,"LRT"=c(object@LRT,NA)))
	# cat("Score test:\n")
	# print(object@scoreT)
	# cat("Likelihood Ratio Test:\n")
	# print(object@LRT)
	})

### more extensive summary...
setMethod("summary",signature("isTest"),function(object){
	
	})



setMethod("plot",signature("isTest"),function(x,y,...){
	par(mfrow=c(1,2))
	plot(x@logLikUnRestr,xlab=expression(beta/xi),ylab=expression(log(f(beta,xi))),...)
	plot(x@logLikRestr,xlab=expression(beta),ylab=expression(log(f(beta,0.5))),...)
	})