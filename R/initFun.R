amisInit <- function(maxit=5000,maxVar=100,s=sqrt(maxVar)){
	
	function(N0,p,target,proposal,verbose=FALSE){
		
		u <- matrix(runif(N0*p),N0,p)
		lu <- log(u/(1-u))
		
		prop0 <- exp(rowSums(lu-2*log(1+exp(lu))))
		
		
		ess <- function(lu,target,prop0){
			function(q){
				N0 <- nrow(lu)
				p <- ncol(lu)
				varini <- s*exp(q)/(1+exp(q))
				xx <- matrix(rep(as.numeric(varini),N0),N0,p,byrow=TRUE)*lu
				targ <- target(xx)
				prop <- prop0/prod(varini)
				w <- exp(targ)/prop
				w[(exp(targ)<sqrt(.Machine$double.eps))&(prop<sqrt(.Machine$double.eps))] <- 0
				-N0/(1+var(w)/(mean(w))^2)
				}
			}

		if(p<=2){
			depo <- c(sqrt(maxVar),rep(1,max(0,p-1)))*sqrt(3/pi^2)
			}else{
				depo <- c(sqrt(maxVar),sqrt(p-2),rep(1,max(0,p-2)))*sqrt(3/pi^2)
				print(depo)
				}
		opto <- optim(log(depo/(s-depo)),ess(lu,target,prop0),control=list(maxit=maxit))
		
		
		#if(verbose) print(opto$counts)
		
		varo <- opto$par
		varini <- s*exp(varo)/(1+exp(varo))
		#if(verbose)print(varini^2*pi^2/3)
		
		xx <- matrix(rep(as.numeric(varini),N0),N0,p,byrow=TRUE)*lu
		targ <- target(xx)
		prop <- prop0/prod(varini)
		w <- exp(targ)/prop
		w[(exp(targ)<sqrt(.Machine$double.eps))&(prop<sqrt(.Machine$double.eps))] <- 0
		
		list(w=w,xx=xx,targ=targ,prop=prop)
			
		}
		
	}


varInit <- function(Var){
	
	function(N0,p,target,proposal,verbose=FALSE){
		
		xx <- proposal$r(N0,mu=rep(0,p),Sig=diag(rep(1,p))) #standardized sample from proposal
		
		ess1 <- unlist(lapply(1:length(Var),function(i){
			xxnew <- sqrt(Var[i])*xx
			targ <- target(xxnew) # returns log
			prop <- proposal$d(xxnew,mu=rep(0,p),Sig=diag(Var[i],p))
			w <- exp(targ)/prop # ratio of densities
			w[(exp(targ)<sqrt(.Machine$double.eps))&(prop<sqrt(.Machine$double.eps))] <- 0
			
			N0/(1+var(w)/(mean(w))^2) 
			}))
		
		vare <- Var[order(ess1,decreasing=TRUE)[1]] #picks which.max(ess1)
		
		# if(verbose){
			# print("vare")
			# print(vare)
			# }
			
		xx <- sqrt(vare)*xx
		targ <- target(xx)
		prop <- proposal$d(xx,mu=rep(0,p),Sig=diag(vare,p))
		w <- exp(targ)/prop
		
		list(w=w,xx=xx,var=vare,targ=targ,prop=prop)
		
		}
		
	}
	
	
uniInit <- function(){
	
	function(N0,p,target,proposal,verbose=FALSE){
		
		u <- matrix(runif(N0*p),N0, p)
		xx <- log(u/(1-u))
		prop <- exp(rowSums(log(u*(1-u))))
		targ <- target(xx)
		w <- exp(targ)/prop
		w[(exp(targ)<sqrt(.Machine$double.eps))&(prop<sqrt(.Machine$double.eps))] <- 0
		
		list(w=w,xx=xx,var=1,targ=targ,prop=prop)
		
		}
	
	}
