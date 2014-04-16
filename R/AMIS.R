AMIS <- function(N,niter,p,target,proposal=mvtComp(df=3),initialize=uniInit(),mixture,verbose=FALSE,parallel = c("no", "multicore", "snow"),nCores=-1,cl=NULL,tol=0.001,seed=NULL,...){
	
	if(is.null(seed)){
		runif(1)
		seed<-sample(.Random.seed,1)
		}
	
	### Get the call and store the argument the function has been called with
	call <- match.call()
	argz <- lapply(2:length(call),function(i)eval(call[[i]]))
	names(argz) <- names(call)[2:length(names(call))]
	
	### parallelization
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    
    if (parallel != "no" && nCores!=1) {
    	require("parallel")
    	
    	if(nCores<1){
    		nCores <- detectCores()
    		
    		}
    	
        if (parallel == "multicore"){
        	have_mc <- .Platform$OS.type != "windows"
        	#if(have_mc)cl <- makeCluster(getOption("cl.cores", nCores))
            }else if (parallel == "snow") have_snow <- TRUE
        
        if (!have_mc && !have_snow) nCores <- 1L
        
        }else nCores <- 1L
	
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
	
	### initial Importance Sample
	J <- sample(1:N[1],N[1],prob=init$w,replace=TRUE)
		
	### from and to determine which row we are adding with the current iteration
	from <- 1
	to <- N[1]
	
	### total Sample Size
	totSS <- N[3]*niter[2]+N[2]*niter[1]+N[1] 
	
	### Matrix with the Importance sample
	xx <- matrix(0,totSS,p)
	xxnew <- init$xx
	xx[from:to,] <- xxnew

	### target, proposal and weights
	w <- prop <- targ <- rep(0, totSS)
	targ[from:to] <- init$targ
	prop[from:to] <- init$prop
	
	### compute the weights
	w[1:to] <- exp(targ[1:to])/prop[1:to]
	### against underflow
	w[1:to][(exp(targ[1:to])<sqrt(.Machine$double.eps))&(prop[1:to]<sqrt(.Machine$double.eps))] <- 0
	
	prop[from:to] <- prop[from:to] * N[1]
	
	### Perplexity
	Perp <- rep(0,sum(niter))
	wBar <- w[from:to]/sum(w[from:to])
	wBar <- wBar[wBar>=sqrt(.Machine$double.eps)] 
	Perp[1] <- exp(-sum(wBar*log(wBar)))/length(wBar)
	
	if(verbose) cat("\nInitialization: local ESS =",N[1]/(1+var(w[1:to])/(mean(w[1:to]))^2),"\n Perplexity =",Perp[1],"\n\n") 
	
	
	### compute the IS estimates of mean and variance
	calcu <- cov.wt(xxnew,wt=w[1:to])
	
	### just one component in the first phase
	alpha <- list()
	alpha[[1]] <- 1
	muHat <- t(calcu$center)
	SigmaHat <- array(calcu$cov,c(dim(calcu$cov),1))
		
	muHat <- list(muHat)
	SigmaHat <- list(SigmaHat)
	
	from <- to + 1
	to <- to + N[2]
	
	cond <- FALSE
	t <- 1
	
	### Adaptation Phase
	if(niter[1]!=0){
		
		if (nCores > 1L && (have_mc || have_snow)) {
			
			sampleImp <- function(n){
				
				### new sample from the proposal
				x <- rprop(n,mu=muHat[[t]][1,],Sig=SigmaHat[[t]][,,1])
				
				### evaluate the target at xxnew
				trgt <- target(x)
				
				### evaluate the proposal at xxnew
				prp <- rep(0,nrow(x))
				prp <- prp + Reduce("+",lapply(1:t,function(k){
					dprop(x,muHat[[k]][1,], SigmaHat[[k]][,,1])
					}))
				ww <- exp(trgt)/prp
				ww[(exp(trgt)<sqrt(.Machine$double.eps))&(exp(prp)<sqrt(.Machine$double.eps))] <- 0
				
				return(list(x=x,trgt=trgt,prp=prp,ww=ww))
				
				}
			
			assign("target",target,environment(sampleImp))
			assign("dprop",dprop,environment(sampleImp))
			assign("rprop",rprop,environment(sampleImp))
							
			}

		
		
		
		
		while((t <=niter[1])&!cond){
			
			
			### Update the proposal
			prop[1:(from-1)] <- prop[1:(from-1)] + N[2]*dprop(xx[1:(from-1),],muHat[[t]][1,],SigmaHat[[t]][,,1])
			### Update the weights
			w[1:(from-1)] <- exp(targ[1:(from-1)])/prop[1:(from-1)]
			### against underflow
			w[1:(from-1)][(exp(targ[1:(from-1)])<sqrt(.Machine$double.eps))&(prop[1:(from-1)]<sqrt(.Machine$double.eps))] <- 0
			
			if (nCores > 1L && (have_mc || have_snow)) {
				
				assign("muHat",muHat,environment(sampleImp))
				assign("SigmaHat",SigmaHat,environment(sampleImp))
				
				if (have_mc){
					
					if (is.null(cl)) {
						cl <- makeCluster(getOption("cl.cores", nCores))
						if (RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
						}
						
					#new <- parallel::mclapply(rep(N[2]/nCores,nCores),sampleImp, ..., mc.cores = nCores)
					jobs <- lapply(rep(ceiling(N[2]/nCores),nCores), function(n) parallel::mcparallel(sampleImp(n), name = n))
					new <- mccollect(jobs)
					#new <- clusterCall(cl,function(n)sampleImp(n),rep(N[2]/nCores,nCores))
					
					
					}else if (have_snow) {
						
						if (is.null(cl)) {
							cl <- parallel::makePSOCKcluster(rep("localhost", nCores))#PSOCK
							#cl <- parallel::makePSOCKcluster(getOption("cl.cores", nCores))
							#cl <- makeCluster(getOption("cl.cores", nCores))
							if (RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
							
							}
						new <- parallel::parLapply(cl,rep(ceiling(N[2]/nCores),nCores),sampleImp)
						
						}
				
				xx[from:to,] <- do.call(rbind,lapply(new,function(n)n$x))[1:length(from:to),]
				targ[from:to] <- unlist(lapply(new,function(n)n$trgt))[1:length(from:to)]
				prop[from:to] <- unlist(lapply(new,function(n)n$prp))[1:length(from:to)]
				w[from:to] <- unlist(lapply(new,function(n)n$ww))[1:length(from:to)]		
					
				}else{
					
					### new sample from the proposal
					xxnew <- rprop(N[2],mu=muHat[[t]][1,],Sig=SigmaHat[[t]][,,1])
					xx[from:to,] <- xxnew  
					
					### evaluate the target at xxnew
					targ[from:to] <- target(xxnew)
					
					### evaluate the proposal at xxnew
					propnew <- rep(0,nrow(xxnew))
					propnew <- propnew + Reduce("+",lapply(1:t,function(k){
						dprop(xxnew,muHat[[k]][1,], SigmaHat[[k]][,,1])
						}))
					prop[from:to] <- propnew
					
					w[from:to] <- exp(targ[from:to])/prop[from:to]
					w[from:to][(exp(targ[from:to])<sqrt(.Machine$double.eps))&(exp(prop[from:to])<sqrt(.Machine$double.eps))] <- 0
					
					}
					
			
			### Compute new IS estimates of Mean and Variance
			calcu <- cov.wt(xx[1:to,],wt=w[1:to])
			SigmaHat[[t+1]] <- array(calcu$cov,c(dim(calcu$cov),1))
			muHat[[t+1]] <- t(calcu$center)
			alpha[[t+1]] <- 1
			
			
			### Compute the perplexity
			wBar <- w[from:to]/sum(w[from:to])
			wBar <- wBar[wBar>=sqrt(.Machine$double.eps)] 
			Perp[t+1] <- exp(-sum(wBar*log(wBar)))/length(wBar)
			
			### Check for convergence
			#cond <- (abs(Perp[t+1]-Perp[t])<=tol[1]*(abs(Perp[t])+tol[2]))
			cond <- (abs(Perp[t+1]-Perp[t])<=tol)
			#cond <- (sqrt(sum((muHat[[t+1]]-muHat[[t]])^2))<=tol)
			
			if(verbose) cat("\nPhase 1: t =",t,"\n Local ESS =",N[2]/(1+var(w[1:to])/(mean(w[1:to]))^2),"\n Perplexity =",Perp[t+1],"\n")
			
			
			from <- to + 1
			to <- to + N[2]
			
			t <- t + 1L
			
			}
		
		}
	
	### Second Phase
	if(niter[2]!=0){
		
		### Draw a sample from the importance distribution
		J <- sample(1:(from-1),N[3],prob=w[1:(from-1)],replace=TRUE)
		
		### Rao-Blackwellized Clustering
		clustMix <- mixture(xx[J,])
		
		### Components of the mixture
		alpha[[niter[1]+1]] <- clustMix$alpha
		muHat[[niter[1]+1]] <- clustMix$muHat
		SigmaHat[[niter[1]+1]] <- clustMix$SigmaHat
		
		to <- to + (N[3]-N[2])
		niter[1] <- min(niter[1],t)

		cond <- FALSE
		niter[1] <- min(niter[1],t)
		
		if (nCores > 1L && (have_mc || have_snow)) {
			
			sampleImp <- function(n){
				
				### Sample xxnew from the mixture...
				x <- do.call(rbind,lapply(1:n,function(j){
					compo <- sample(1:G,1,prob=alpha[[t]])
					y <- t(t(rprop(1,muHat[[t]][compo,],SigmaHat[[t]][,,compo])))
					rownames(y) <- as.character(compo)
					y
					}))
				
				### evaluate the proposal at xxnew
				prp <- rep(0,nrow(x))
				prp <- prp + Reduce("+",lapply(1:t,function(k){Reduce("+",lapply(1:length(alpha[[k]]),function(g){
						(alpha[[k]][g]* dprop(x,muHat[[k]][g,],SigmaHat[[k]][,,g]))
						}))}))
					
				### evaluate the target at xxnew
				trgt <- target(x)
				ww <- exp(trgt)/(prp)
				ww[(exp(trgt)<sqrt(.Machine$double.eps))&(prp<sqrt(.Machine$double.eps))] <- 0
				
				return(list(x=x,trgt=trgt,prp=prp,ww=ww))
				
				}
				
				
			assign("target",target,environment(sampleImp))
			assign("dprop",dprop,environment(sampleImp))
			assign("rprop",rprop,environment(sampleImp))
				
			}
		
		while((t <= niter[2]+niter[1])&!cond){
			
			newEnv <-  new.env()
			
			G <- length(alpha[[t]])
			
			
			### Update the proposal
			prop[1:(from-1)] <- prop[1:(from-1)] + Reduce("+",lapply(1:G,function(g){
				N[3]*(alpha[[t]][g]* dprop(xx[1:(from-1),],muHat[[t]][g,],SigmaHat[[t]][,,g]))
				}))
			### Update the weights
			w[1:(from-1)] <- exp(targ[1:(from-1)])/prop[1:(from-1)]
			### against underflow
			w[1:(from-1)][(exp(targ[1:(from-1)])<sqrt(.Machine$double.eps))&(prop[1:(from-1)]<sqrt(.Machine$double.eps))] <- 0
			
			if (nCores > 1L && (have_mc || have_snow)) {
				
				assign("alpha",alpha,environment(sampleImp))
				assign("muHat",muHat,environment(sampleImp))
				assign("SigmaHat",SigmaHat,environment(sampleImp))
				
				if (have_mc){
					
					if (is.null(cl)) {
						cl <- makeCluster(getOption("cl.cores", nCores))
						if (RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
						}
						
					
					jobs <- lapply(rep(ceiling(N[3]/nCores),nCores), function(n) mcparallel(sampleImp(n), name = n))
					new <- mccollect(jobs)
					
					}else if (have_snow) {
						#list(...)
						if (is.null(cl)) {
							cl <- parallel::makePSOCKcluster(rep("localhost", nCores))
							#cl <- parallel::makePSOCKcluster(nCores)
							#cl <- parallel::makePSOCKcluster(getOption("cl.cores", nCores))
							if (RNGkind()[1L] == "L'Ecuyer-CMRG") parallel::clusterSetRNGStream(cl)
							#parallel::clusterEvalQ(cl, library(survival))
							new <- parallel::parLapply(cl,rep(ceiling(N[3]/nCores),nCores),sampleImp)
							#parallel::stopCluster(cl)
							
							}else{
								#parallel::clusterEvalQ(cl, library(survival))
								new <- parallel::parLapply(cl,rep(N[3]/nCores,nCores),sampleImp)
								}
							}
							
				xx[from:to,] <- do.call(rbind,lapply(new,function(n)n$x))[1:length(from:to),]
				targ[from:to] <- unlist(lapply(new,function(n)n$trgt))[1:length(from:to)]
				prop[from:to] <- unlist(lapply(new,function(n)n$prp))[1:length(from:to)]
				w[from:to] <- unlist(lapply(new,function(n)n$ww))[1:length(from:to)]		
				
					
				}else{
					
					
					### Sample xxnew from the mixture...
					xxnew <- do.call(rbind,lapply(1:N[3],function(j){
						compo <- sample(1:G,1,prob=alpha[[t]])
						x <- t(t(rprop(1,muHat[[t]][compo,],SigmaHat[[t]][,,compo])))
						rownames(x) <- as.character(compo)
						x
						}))
					### Update the proposal
					prop[1:(from-1)] <- prop[1:(from-1)] + Reduce("+",lapply(1:G,function(g){
						N[3]*(alpha[[t]][g]* dprop(xx[1:(from-1),],muHat[[t]][g,],SigmaHat[[t]][,,g]))
						}))
					
					propnew <- rep(0,nrow(xxnew))
					propnew <- propnew + Reduce("+",lapply(1:t,function(k){Reduce("+",lapply(1:length(alpha[[k]]),function(g){
							(alpha[[k]][g]* dprop(xxnew,muHat[[k]][g,],SigmaHat[[k]][,,g]))
							}))}))
					
					xx[from:to,] <- xxnew  
					
					targ[from:to] <- target(xxnew)
					prop[from:to] <- propnew
					
					w[from:to] <- exp(targ[from:to])/(prop[from:to])
					w[from:to][(exp(targ[from:to])<sqrt(.Machine$double.eps))&(prop[from:to]<sqrt(.Machine$double.eps))] <- 0
					
					}
					
			J <- sample(1:to,N[3],prob=w[1:to],replace=TRUE)
			clustMix <- mixture(xx[J,])
			
			alpha[[t+1]] <- clustMix$alpha
			muHat[[t+1]] <- clustMix$muHat
			SigmaHat[[t+1]] <- clustMix$SigmaHat
			
			wBar <- w[from:to]/sum(w[from:to])
			wBar <- wBar[wBar>=sqrt(.Machine$double.eps)] 
			Perp[t+1] <- exp(-sum(wBar*log(wBar)))/length(wBar)
						
			if(verbose)cat("\nPhase 2: t =",t,"Number of Components G =",G,"\n Local ESS =",N[3]/(1+var(w[1:to])/(mean(w[1:to]))^2),"\n Perplexity =",Perp[t+1],"\n")
			
			#cond <- (abs(Perp[t+1]-Perp[t])>=tol[1]*(abs(Perp[t])+tol[2]))
			cond <- (abs(Perp[t+1]-Perp[t])<=tol)
			#cond <- (sqrt(sum((alpha[[t+1]]%*%muHat[[t+1]]-alpha[[t]]%*%muHat[[t]])^2))<=tol)
			
			from <- to+1
			to <- to + N[3]
			
			t <- t + 1L	
			
			}
			
		niter[2] <- min(niter[2],t)
		
		}
		
		
	### Delete the Clusters
	if(have_snow||have_mc)parallel::stopCluster(cl)
	
	### Select the mixture with the highest perplexity	
	best <- which.max(Perp)
	maxPerp <- max(Perp)
	
	IS <- cbind(xx,w)[1:(from-1),]
	colnames(IS) <- c(paste("x",1:p,sep=""),"w")
	
	ESS <- nrow(IS)/(1+var(IS[,p+1])/(mean(IS[,p+1]))^2)
	
	G <- length(alpha[[best]])
	
	if(verbose)cat("\nGlobal ESS =",ESS,"\nMaximum Perplexity =",maxPerp,"\n",G,"Components\n\n")
	
	
	### this environment will contain the d-, p- and r- functions relative to the target.
	env <- new.env()
	
	### the prefix has the same meaning as in dnorm, pnorm and rnorm.
	dTarg <- function(pp=alpha[[best]],mu=muHat[[best]],Sig=SigmaHat[[best]]){
		G <- length(pp)
		function(xx,log=TRUE){
			
			px <- sum(exp(unlist(lapply(1:G,function(g)log(pp[g])+dprop(xx,mu[g,],Sig[,,g],log=TRUE)))))
			if(log) px <- log(px)
			
			px
			
			}
		}
	dTarg <- dTarg()
	
	rTarg <- function(IS=IS,pp=alpha[[best]],mu=muHat[[best]],Sig=SigmaHat[[best]]){
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
	if(name=="") name <- "Targ"
	
	assign(paste("d",name,sep=""),dTarg,envir=env)
	assign(paste("r",name,sep=""),rTarg,envir=env)
	assign(paste("p",name,sep=""),pTarg,envir=env)
	
	argz[["niter"]] <- niter
	
	return(new("ISO",IS=IS,
	                 Prop=alpha[[best]],
	                 Mean= muHat[[best]],
	                 Var=SigmaHat[[best]],
	                 ESS=ESS,
	                 Perp = maxPerp,
	                 seed=seed,
	                 envir=env,
	                 call=call,
	                 args=argz))	
	}

#amis <- AMIS(N=1000,niter=10,p=2,target=targetBanana(b=0.1),initialize= amisInit(maxit=5000),mixture=mclustMix(),verbose=TRUE)	
