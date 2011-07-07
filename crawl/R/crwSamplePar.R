crwSamplePar <- function(object.sim, method="IS", size=1000, df=Inf, grid.eps=1, crit=2.5, scale=1)
{
   if(!inherits(object.sim, 'crwSimulator'))
      stop("Argument needs to be of class 'crwSimulator'\nUse 'crwSimulator( )' to create")
   fixPar <- object.sim$fixPar
   Cmat <- object.sim$Cmat[is.na(fixPar),is.na(fixPar)]
   se <- sqrt(diag(Cmat))
   err.mfX <- object.sim$err.mfX
   err.mfY <- object.sim$err.mfY
   parMLE <- object.sim$par
   n2ll.mode <- -2*object.sim$loglik
   stopMod <- object.sim$stopMod
   driftMod <- object.sim$driftMod
   stop.mf <- object.sim$stop.mf
   mov.mf <- object.sim$mov.mf
   y <- object.sim$y
   x <- object.sim$x
   loctype <- object.sim$loctype
   delta <- object.sim$delta
   a1.y <- object.sim$a1.y
   a1.x <- object.sim$a1.x
   P1.x <- object.sim$P1.x
   P1.y <- object.sim$P1.y
   lonAdj <- object.sim$lonAd
   n.errX <- object.sim$n.errX
   n.errY <- object.sim$n.errY
   n.mov <- object.sim$n.mov
   N <- object.sim$N
   lower <- object.sim$lower
   upper <- object.sim$upper 
   prior <- object.sim$prior
   cat("\nComputing importance weights ...\n")
   if(method=="IS"){
       thetaMat <- matrix(NA, size, length(fixPar)+3)
	   for(i in 1:(size-1)){
   	 	par <- parMLE
     	eInd <- is.na(fixPar)
   	 	eps <- rmvtt(mu=rep(0,sum(eInd)), Sigma=scale*Cmat, df=df, lower-par[eInd], upper-par[eInd])
     	par[eInd] <- parMLE[eInd] + eps
     	if(df==Inf) dens <- dmvnorm(eps, sigma=scale*Cmat, log=TRUE) - dmvnorm(0.0*eps, sigma=scale*Cmat, log=TRUE)
     	else dens <- dmvt(eps, sigma=scale*Cmat, df=df, log=TRUE) - dmvt(0.0*eps, sigma=scale*Cmat, df=df, log=TRUE)
		n2ll.val <- crwN2ll(par[eInd], fixPar, y, x, loctype, delta, a1.y, a1.x,
                      					P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
                      					n.errX, n.errY, n.mov, stopMod, driftMod, prior=prior, need.hess=FALSE)
		thetaMat[i,] <- c(-n2ll.val/2 - dens, -n2ll.val/2, dens, par)
	  }
	  thetaMat[size,] <- c(object.sim$loglik, object.sim$loglik, 0, object.sim$par)
	  thetaMat[,1] <- exp(thetaMat[,1]-max(thetaMat[,1]))/sum(exp(thetaMat[,1]-max(thetaMat[,1])))
	}
else if(method=="mcmc"){

}
else if(method=="quadrature"){
	Eigen.list <- eigen(Cmat, symmetric=TRUE)
	V <- Eigen.list$vectors
	D <- diag(sqrt(Eigen.list$values))
	np <- sum(is.na(fixPar))
	grid.list <- rep(list(0), np)
	eInd <- is.na(fixPar)
	thetaMat <- matrix(c(-n2ll.mode/2, -n2ll.mode/2, 0, parMLE), nrow=1)
	   for(k in 1:np){
	   		stop.grid <- TRUE
	   		z <- rep(0,np)
	   		while(stop.grid){
	   			z[k] <- z[k] + grid.eps 
		 		par <- parMLE
		 		par[eInd] <- parMLE[eInd] + V%*%D%*%z
		 		if(any(par[eInd]>upper) | any(par[eInd]<lower)) stop.grid <- FALSE
		 		else{
		 		 	n2ll.val <- crwN2ll(par[eInd], fixPar, y, x, loctype, delta, a1.y, a1.x,
                      					P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
                      					n.errX, n.errY, n.mov, stopMod, driftMod, prior=prior, need.hess=FALSE)
		 		 	if(-(n2ll.mode - n2ll.val)/2 > crit) stop.grid <- FALSE
		 		 	else{
		 		 		grid.list[[k]] <- c(grid.list[[k]],z[k])
		 		 		thetaMat <- rbind(thetaMat, c(-n2ll.val/2, -n2ll.val/2, 0, par))
		 		 	}
		 		 }
		 	}
	   		stop.grid <- TRUE
	   		z <- rep(0,np)
	   		while(stop.grid){
	   			z[k] <- z[k] - grid.eps 
		 		par <- parMLE
		 		par[eInd] <- parMLE[eInd] + V%*%D%*%z
		 		if(any(par[eInd]>upper) | any(par[eInd]<lower)) stop.grid <- FALSE
		 		else{
		 		 	n2ll.val <- crwN2ll(par[eInd], fixPar, y, x, loctype, delta, a1.y, a1.x,
                      					P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
                      					n.errX, n.errY, n.mov, stopMod, driftMod, prior=prior, need.hess=FALSE)
		 		 	if(-(n2ll.mode - n2ll.val)/2 > crit) stop.grid <- FALSE
		 		 	else{
		 		 		grid.list[[k]] <- c(grid.list[[k]],z[k])
		 		 		thetaMat <- rbind(thetaMat, c(-n2ll.val/2, -n2ll.val/2, 0, par))
		 		 	}
		 		 }
		 	}
	  }
	  grid.pts <- as.matrix(expand.grid(grid.list))
	  grid.pts <- grid.pts[apply(grid.pts==0, 1, sum) < np-1, ]
	  numEvals <- nrow(grid.pts)+nrow(thetaMat)
	  cat("\nEvaluating ", nrow(grid.pts)+nrow(thetaMat), " quadrature points ...\n")
	  parFix <- ifelse(!eInd, parMLE, 0)
#	  grid.eval <- system.time(t(apply(grid.pts, 1, 
#			  			function(z,...){
#							par[eInd] <- parMLE[eInd] + V%*%D%*%z
#							if(any(par[eInd]>upper) | any(par[eInd]<lower)) return(rep(NA,length(par)+3))
#							else{
#								n2ll.val <- crwN2ll(par[eInd], fixPar, y, x, loctype, delta, a1.y, a1.x,
#										P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
#										n.errX, n.errY, n.mov, stopMod, driftMod)
#								if(-(n2ll.mode - n2ll.val)/2 > crit) return(rep(NA,length(par)+3))
#								else c(-n2ll.val/2, -n2ll.val/2, 0, par)
#							}
#						},fixPar, y, x, loctype, delta, a1.y, a1.x,
#						P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
#						n.errX, n.errY, n.mov, stopMod, driftMod)))
	 for(i in 1:nrow(grid.pts)){
	  		   	z <- grid.pts[i,]
		 		par <- parMLE
		 		par[eInd] <- parMLE[eInd] + V%*%D%*%z
		 		if(any(par[eInd]>upper) | any(par[eInd]<lower)) next
		 		else{
		 		 	n2ll.val <- crwN2ll(par[eInd], fixPar, y, x, loctype, delta, a1.y, a1.x,
                      					P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
                      					n.errX, n.errY, n.mov, stopMod, driftMod, prior=prior, need.hess=FALSE)
		 		 	if(-(n2ll.mode - n2ll.val)/2 > crit) next
		 		 	else thetaMat <- rbind(thetaMat, c(-n2ll.val/2, -n2ll.val/2, 0, par))
		 		 }
			 }
thetaMat[,1] <- exp(thetaMat[,1]-max(thetaMat[,1]))/sum(exp(thetaMat[,1]-max(thetaMat[,1])))
}
else stop("\nIncorrect specification of parameter sampling method\n")

colnames(thetaMat) <- c("w", "lik", "prop.lik", object.sim$nms)
attr(thetaMat,"effSamp") <- nrow(thetaMat)/(1+(sd(thetaMat[,"w"])/mean(thetaMat[,"w"]))^2) 
attr(thetaMat, "method") <- method
attr(thetaMat, "numLikEval") <- ifelse(method=="quadrature", numEvals, size)
if(is.null(object.sim$thetaSampList)) object.sim$thetaSampList <- list(thetaMat)
else object.sim$thetaSampList <- append(object.sim$thetaSampList, list(thetaMat))
return(object.sim)	
}
     