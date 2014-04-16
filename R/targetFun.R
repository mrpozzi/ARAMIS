### Banana Shaped Target
## Moderate banana shape: b=0.01 and sig=100
## Strong banana shape: b=0.03 and sig=100

targetBanana <- function(sig=100,b=0.01,log=TRUE){
	
	
	function(x){
		if (is.vector(x)) x <- t(as.matrix(x))
		p <- ncol(x)
		x[,2] <- x[,2]+b*x[,1]^2-b*sig
		
		Mu <- rep(0,p)
		Sigma <- diag(p)
		Sigma[1,1] <- sig
		
		dmnorm(x,mean=Mu,varcov=Sigma,log=log)
		
		}
	
	}
	
### Multivariate Normal Target

targetMVN <- function(Mu=rep(0,2),Sigma=diag(1,2,2)){
	
	function(x){
		p <- length(Mu)
		-p/2*(log(2*pi))-1/2*log(det(Sigma))-1/2*t(x-Mu)%*%ginv(Sigma)%*%(x-Mu)
		}
		
	}



### Mixture Target

targetMix <- function(alpha=rep(.5,2),Mu=matrix(0,2,2),Sigma=array(cbind(diag(1,2,2),diag(1,2,2)),dim=c(2,2,2))){
	function(x){
		if(sum(alpha)!=1)alpha <- alpha/sum(alpha)
		Reduce("+",lapply(1:length(alpha),function(g){
			if(length(dim(x))>1){
				alpha[g]*apply(x,1,function(y)dmnorm(y,Mu[g,],Sigma[,,g]))
				}else{
					log(alpha[g]*dmnorm(x,Mu[g,],Sigma[,,g]))
					}
			}))
		}
	}

