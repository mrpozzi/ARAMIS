mvtComp <- function(df=3){
	list("d"=function(xx,mu=rep(0,ncol(xx)),Sig=diag(1,ncol(xx),ncol(xx)),log=FALSE){
		dmt(xx,mean=mu,S=Sig,df=df,log=log)
		},"r"=function(n=1,mu=0,Sig=1){
			rmt(n,mean=mu,S=Sig,df=df)
			})
	}
	
mvnormComp <- function(){
	list("d"=function(xx,mu=rep(0,ncol(xx)),Sig=diag(1,ncol(xx),ncol(xx)),log=FALSE){
		dmnorm(xx,mean=mu,varcov=Sig,log=log)
		},"r"=function(n=1,mu=0,Sig=1){
		rmnorm(n,mean=mu,varcov=Sig)
		})
	}

