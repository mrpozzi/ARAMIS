AIS <- function(N,niter,p,target,proposal=mvtComp(df=3),initialize=uniInit(),mixture,verbose=FALSE,tol=0.001,seed=NULL){
	
	if(is.null(seed)){
		runif(1)
		seed<-sample(.Random.seed,1)
		}
	
	### Get the call and store the argument the function has been called with
	call <- match.call()
	argz <- lapply(2:length(call),function(i)eval(call[[i]]))
	names(argz) <- names(call)[2:length(names(call))]
	
	### make sure the number of groups is acceptable
	G <- get("G",envir=environment(mixture))
	if(any(G > min(N))) G <- G[G <= min(N)]
	assign("G",sort(G),envir=environment(mixture))
	
	### set the RNG seed
	set.seed(seed)
	
	### same Sample Size
	if(length(N)==1) N <- rep(N,3)
	### same number of iterations
	if(length(niter)==1) niter <- rep(niter,2)
	
	### Get the density and the sampling fuction for the component of the proposal
	dprop <- proposal$d
	rprop <- proposal$r
	
	### Set the initial condition
	init <- initialize(N[1],p,target,proposal,verbose)
	
	### from and to determine which row we are adding with the current iteration
	from <- 1
	to <- N[1]
	
	### total Sample Size
	totIt <- N[3]*niter[2]+N[2]*niter[1]+N[1]
	
	### Matrix with the Importance sample and the weights
	IS <- matrix(0,nrow= totIt,ncol=p+1)
	IS[1:N[1],] <- cbind(init$xx,init$w)
	
	### Perplexity
	Perp <- rep(0,sum(niter))
	wBar <- init$w/sum(init$w)
	wBar <- wBar[wBar>=sqrt(.Machine$double.eps)]
	maxPerp <- Perp[1] <- exp(-sum(wBar*log(wBar)))/length(wBar)
	
	if(verbose)cat("\nInitialization: local ESS =",N[1]/(1+var(IS[,p+1])/(mean(IS[,p+1]))^2),"\n Perplexity =",Perp[1],"\n\n") 
	
	from <- to+1
	to <- to + N[2]
	
	t <- 1L
	cond <- FALSE
	
	### Adaptation Phase
	if(niter[1]!=0){
	
		while((t <= niter[1])&!cond){
			
			### Compute new IS estimates of Mean and Variance
			calcu <- cov.wt(IS[1:(from-1),1:p],wt=IS[1:(from-1),p+1])
			Sigma <- calcu$cov
			Mean <- calcu$center
			
			### new sample from the proposal
			xx <- rprop(N[2],mu=Mean,Sig=Sigma)
			
			### evaluate the target at xx
			targ <- target(xx)
			
			### evaluate the proposal at xxnew
			prop <- dprop(xx,mu= Mean,Sig=Sigma)
			w <- exp(targ)/prop
			w[(exp(targ)<sqrt(.Machine$double.eps))&(prop<sqrt(.Machine$double.eps))] <- 0
			
			### Compute the perplexity
			wBar <- w/sum(w)
			wBar <- wBar[wBar>=sqrt(.Machine$double.eps)] 
			Perp[t+1] <- exp(-sum(wBar*log(wBar)))/length(wBar)
	
			
			if(verbose) cat("\nPhase 1: t =",t,"\n Local ESS =",N[2]/(1+var(w)/(mean(w))^2),"\n Perplexity =",Perp[t+1],"\n") 
				
	
			IS[from:to,1:ncol(IS)] <- cbind(xx,w)
			from <- to + 1
			to <- to + N[2]
			
			### Check for convergence
			#cond <- (abs(Perp[t+1]-Perp[t])<=tol[1]*(abs(Perp[t])+tol[2]))
			cond <- (abs(Perp[t+1]-Perp[t])<=tol)
			#cond <- (sqrt(sum((muHat[[t+1]]-muHat[[t]])^2))<=tol)
			
			t <- t + 1L
			
			if(Perp[t+1]>maxPerp){
				maxPerp <- Perp[t+1]
				muHatB <- Mean
				SigmaHatB <- Sigma
				alphaB <- 1
				}
			
			}
		
		if(verbose)cat("\nGlobal ESS =",N[2]/(1+var(w)/(mean(w))^2),"\nPerplexity =",Perp[t],"\n") 
		
		}
		
	### Draw a sample from the importance distribution
	J <- sample(1:(from-1),N[2],prob=IS[1:(from-1),p+1],replace=TRUE) 
	xx <- IS[J,1:p]
	
	to <- to + (N[3]-N[2]) 
	
	### Second Phase
	if (niter[2]!=0){
		
		t <- 1
		cond <- FALSE
		
		while((t <= niter[2])&!cond){
			
			### Rao-Blackwellized Clustering
			clustMix <- mixture(xx)
			G <- clustMix$G
			cluster <- clustMix$cluster
			
			### Components of the mixture
			pp <- clustMix$alpha
			muHat <- clustMix$muHat
			varHat <- clustMix$SigmaHat
			
			### Sample xxnew from the mixture...
			xx <- do.call(rbind,lapply(1:N[3],function(j){
				compo <- sample(1:G,1,prob=pp) 
				x <- t(t(rprop(1,muHat[compo,], varHat[,,compo])))
				rownames(x) <- as.character(compo)
				x
				}))
			
			prop <- rep(0,N[3]) # density
			prop <- prop + Reduce("+",lapply(1:G,function(g){
				pp[g]* dprop(xx,mu=muHat[g,],Sig=varHat[,,g])
				}))
			
			targ <- target(xx) # log(density)
			
			w <- exp(targ)/prop
			w[(exp(targ)<sqrt(.Machine$double.eps))&(prop<sqrt(.Machine$double.eps))] <- 0
			
			IS[from:to,] <- cbind(xx,w)
			from <- to + 1
			to <- to + N[3] 
			
			
			J <- sample(1:(from-1),N[3],prob=IS[1:(from-1),p+1],replace=TRUE)
			xx <- IS[J,1:p]
			
			calcu <- cov.wt(IS[1:(from-1),1:p],wt=IS[1:(from-1),p+1])
			Sigma <- calcu$cov
			Mean <- calcu$center
			
			wBar <- w/sum(w)
			wBar <- wBar[wBar>=sqrt(.Machine$double.eps)] 
			Perp[t+1] <- exp(-sum(wBar*log(wBar)))/length(wBar)
			
			
			if(verbose)cat("Phase 2: t =",t,"\nNumber of Components G =",G,"\n Local ESS =",N[3]/(1+var(w)/(mean(w))^2),"\n Perplexity =",Perp[t+1],"\n")
			
			### Select the mixture with the highest perplexity	
			if(Perp[t+1]>maxPerp){
				maxPerp <- Perp[t+1]
				muHatB <- muHat
				SigmaHatB <- varHat
				alphaB <- pp
				}
			
			#cond <- (abs(Perp[t+1]-Perp[t])>=tol[1]*(abs(Perp[t])+tol[2]))
			cond <- (abs(Perp[t+1]-Perp[t])<=tol)
			#cond <- (sqrt(sum((alpha[[t+1]]%*%muHat[[t+1]]-alpha[[t]]%*%muHat[[t]])^2))<=tol)
			
			t <- t + 1L
			
			}
		}
	
	gc()
	
	
	colnames(IS) <- c(paste("x",1:p,sep=""),"w")
	
	ESS <- nrow(IS)/(1+var(IS[,p+1])/(mean(IS[,p+1]))^2)
	
	if(verbose)cat("\nGlobal ESS =",ESS,"\nMaximum Perplexity =",maxPerp,"\n\n")
	
	### this environment will contain the d-, p- and r- functions relative to the target.
	env <- new.env()
	
	### the prefix has the same meaning as in dnorm, pnorm and rnorm.
	dTarg <- function(pp=alphaB,mu=muHatB,Sig=SigmaHatB){
		G <- length(pp)
				function(xx,log=TRUE){
			if(log){
				sum(unlist(lapply(1:G,function(g)log(pp[g])+dprop(xx,mu[g,],Sig[,,g],log=TRUE))))
				}else{
					exp(sum(unlist(lapply(1:G,function(g)log(pp[g])+dprop(xx,mu[g,],Sig[,,g],log=TRUE)))))
					}
			
			}
		}
	dTarg <- dTarg()
	
	rTarg <- function(IS=IS,pp=alphaB,mu=muHatB,Sig=SigmaHatB){
		G <- length(pp)
		function(n,ISmpl=FALSE){
			
			if(ISmpl){
				J <- sample(1:nrow(IS),size=n,prob=IS[,ncol(IS)],replace=TRUE)
				IS[J,-ncol(IS)]
				}else{
					do.call(rbind,lapply(1:n,function(j){
						compo <- sample(1:G,1,prob=pp)
						x <- t(t(rprop(1,mu[compo,], Sig[,,compo])))
						rownames(x) <- as.character(compo)
						x
						}))
					}
			}
		}
	rTarg <- rTarg()
	
	pTarg <- function(IS){
		function(q,lower.tail=TRUE,log.p=FALSE,N=1000){
			J <- sample(1:nrow(IS),N,prob= IS[,ncol(IS)],replace=TRUE)
			p <- sum(apply(IS[J,1:(ncol(IS)-1)],1,function(x)min(x<=q)))/N
			if(!lower.tail)p <- 1-p
			if(log.p) p <- log(p)
			p
			}
		}
	pTarg <- pTarg(IS)
	
	name <- gsub("target","",unlist(strsplit(deparse(call[["target"]]),split="(",fixed=TRUE)[1])[1],ignore.case=TRUE)
	
	assign(paste("d",name,sep=""),dTarg,envir=env)
	assign(paste("r",name,sep=""),rTarg,envir=env)
	assign(paste("p",name,sep=""),pTarg,envir=env)
	
	
	return(new("ISO",IS=IS,
	                 Prop=alphaB,
	                 Mean= muHatB,
	                 Var=SigmaHatB,
	                 ESS=ESS,
	                 Perp = maxPerp,
	                 seed=seed,
	                 envir=env,
	                 call=call,
	                 args=argz))
	
		
	}
	
#ais <- AIS(N=1000,niter=10,p=2,target=targetBanana(b=0.1),initialize= amisInit(maxit=5000),mixture=mclustMix(),verbose=TRUE)	
