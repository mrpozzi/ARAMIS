#################################
### Importance Sampling Object

setClass("ISO",representation=representation(IS = "matrix",
                                             Prop="numeric",
                                             Mean="matrix",
                                             Var = "array",
                                             Perp = "numeric",
                                             ESS="numeric",
                                             seed="integer",
                                             envir="environment",
                                             call="call",
                                             args="list"))


### Minimal description
setMethod("show",signature("ISO"),function(object){
	cat("Call:\n")
	print(object@call)
	})

### more extensive summary...
setMethod("summary",signature("ISO"),function(object){
	cat("Call:\n")
	print(object@call)
	cat("\nGlobal ESS: ")
	cat(object@ESS,"\n")
	cat("Maximum Perplexity: ")
	cat(object@Perp,"\n")
	#cat("Number of Components: ")
	#cat(object@Perp,"\n")
	})

### plot method
setMethod("plot", signature = "ISO", definition = function (x,whichOne=1L:3L,prj=1L:2L,N=100,xlim=c(-30,30),ylim=c(-30,30),main,phi=30, theta=-30,border=NA,xlab=expression(x[1]),ylab=expression(x[2]),zlab="target",n=100,...)  {
	
	show <- rep(FALSE, 5)
    show[whichOne] <- TRUE
	
	if(length(whichOne)>1) par(ask=TRUE)
	if(missing(main)) main <- paste(x@call[[1]],"Algorithm;",x@call[["target"]],"Target",sep=" ")
	
	x1 <- seq(xlim[1],xlim[2],length.out=n)
	x2 <- seq(ylim[1],ylim[2],length.out=n)
	J <- sample(1:nrow(x@IS),N,prob=x@IS[,ncol(x@IS)],replace=TRUE)
	xx <- x@IS[J,prj]
	
	p <- ncol(x@IS)-1
	
	target <- x@args$target
	tag <- matrix(0,n,n)
	y <- rep(0,p)
	for (i in 1:n) for (j in 1:n){
		y[prj] <- c(x1[i],x2[j])
		tag[i,j] <- exp(target(y))
		}
	tag[is.nan(tag)] <- NA
	tag[is.infinite(tag)] <- NA
	
	### image
	if(show[1L]){
		image(x1,x2,tag,col=heat.colors(100),main=main,xlab=xlab,ylab=ylab)#,...)
		points(xx,col="blue",pch=3,...)
		}
			
	###	3-D plot
	if(show[2L]) {
		
    	myColors <- colorRampPalette( c("blue", "green") )
    	nbcol <- 100
    	color <- myColors(nbcol)
    	zfacet <- tag[-1, -1] + tag[-1, -nrow(tag)] + tag[-nrow(tag), -1] + tag[-nrow(tag), -nrow(tag)]
    	facetcol <- cut(zfacet, nbcol)

    	res <- persp(x1,x2,tag, col=color[facetcol], phi=phi,theta=theta,border=border,xlab=xlab,ylab=ylab,zlab=zlab,main=main)
    	
    	
    	clr <- adjustcolor("red",alpha=0.5)
    	points(trans3d(xx[,1], xx[,2], exp(target(xx)), pmat = res), col =clr, pch =20,main=main,xlab=xlab,ylab=ylab)
    	
    	} 
    
    ### contour
    if(show[3L]){    
    	#browser()
    	y <- rep(0,p)
    	dTarg <- get(ls(x@envir)[1],envir=x@envir)
    	den <- matrix(0,n,n)
    	for (i in 1:n) for (j in 1:n){
    		y[prj] <- c(x1[i],x2[j])
    		den[i,j] <- exp(dTarg(y,log=TRUE))
    		}
    		
    	contour(x1,x2,tag, xlim=xlim,ylim=ylim,nlevels=15,col="red")
    	contour(x1,x2,den,nlevels=15,col="blue",add=TRUE)
    	contour(kde2d(xx[,1], xx[,2],lims = c(range(xlim), range(ylim))),nlevels=15,add=TRUE,)
    	}   
    if(show[5L]){
		
		}
							
	
	})
	
### subsetting method
setMethod("[[", signature = "ISO", definition = function (x, i,j,..., drop = TRUE) {
	if (!missing(j)) {stop("Wrong number of dimensions.")}
	if (!missing(i)) {
		return(switch(class(i),"character" = attributes(x)[[i]],
		                       "integer" = attributes(x)[[i]],
		                       "numeric" = attributes(x)[[i]],
		                       "logical" = attributes(x)[c(i,rep(FALSE,length(attributes(x))-length(i)))],
		                        stop("Subsetting object needs to be either a character, a numeric/integer or a logical.")
		                        ))
		}else{return(NULL)}
  
   })

### subsetting method
setMethod("$", signature = "ISO", definition = function(x, name) {
	 x[[name]]
	 })

### alias for slotNames 
setMethod("names", signature = "ISO", definition = function(x) {
	 slotNames(x)
	 })
	 

### returns the mean OR computes the mean of f(X)
setMethod("mean", signature = "ISO", definition = function(x, fun=NULL,...) {
	 if(is.null(fun)){
	 	calcu <- cov.wt(x@IS[,1:(ncol(x@IS)-1)],wt=x@IS[,"w"],...)
	 	}else{
	 		calcu <- cov.wt(fun(x@IS[,1:(ncol(x@IS)-1)]),wt=x@IS[,"w"],...)
	 		}
	 
	 return(calcu$center)
	 
	 })

### returns the variance OR computes the variance of f(X)
setMethod("var", signature = "ISO", definition = function(x, y=NULL,na.rm=FALSE,use) {
	 if(is.null(y)){
	 	calcu <- cov.wt(x@IS[,1:(ncol(x@IS)-1)],wt=x@IS[,"w"])
	 	}else{
	 		calcu <- cov.wt(y(x@IS[,-ncol(x@IS)]),wt=x@IS[,"w"])
	 		}
	 
	 return(calcu$cov)
	 
	 })

### returns the importance weights 
setMethod("weights",signature("ISO"),function(object){
	return(object@IS[,"w"])
	})

## returns an importance sample  
setMethod("sample",signature("ISO"),function(x,size){
	J <- sample(1:nrow(x@IS),size,prob=x@IS[,ncol(x@IS)],replace=TRUE)
	x@IS[J,-ncol(x@IS)]
	})

# setMethod("attach",signature("ISO"),function(what, warn.conflicts = TRUE){#,...
	# attach(what@envir, warn.conflicts= warn.conflicts)#,...
	# })

# setMethod("detach",signature("ISO"),function(name){
	# detach(name@envir)
	# })
