`crwPostIS` <-
function(object.sim, fullPost=TRUE, df=Inf, scale=1, thetaSamp=NULL)
################################################################################
################################################################################
{
   if(!inherits(object.sim, 'crwSimulator'))
      stop("Argument needs to be of class 'crwSimulator'\nUse 'crwSimulator( )' to create")
   fixPar <- object.sim$fixPar
   Cmat <- object.sim$Cmat[is.na(fixPar),is.na(fixPar)]
   se <- sqrt(diag(Cmat))
   err.mfX <- object.sim$err.mfX
   err.mfY <- object.sim$err.mfY
   stopMod <- object.sim$stopMod
   driftMod <- object.sim$driftMod
   stop.mf <- object.sim$stop.mf
   mov.mf <- object.sim$mov.mf
   n.errX <- object.sim$n.errX
   n.errY <- object.sim$n.errY
   n.mov <- object.sim$n.mov
   N <- object.sim$N
   lower <- object.sim$lower
   upper <- object.sim$upper
   loctype <- object.sim$loctype
   par <- object.sim$par
   
   ###
   ### Sample parameter vector
   ###
   if(fullPost & is.null(object.sim$thetaSampList)) {
      eInd <- is.na(fixPar)
      eps <- rmvtt(mu=rep(0,sum(eInd)), Sigma=scale*Cmat, df=df, lower-par[eInd], upper-par[eInd])
      par[eInd] <- par[eInd] + eps
      if(df==Inf) dens <- dmvnorm(eps, sigma=scale*Cmat, log=TRUE) - dmvnorm(0.0*eps, sigma=scale*Cmat, log=TRUE)
      else dens <- dmvt(eps, sigma=scale*Cmat, df=df, log=TRUE) - dmvt(0.0*eps, sigma=scale*Cmat, df=df, log=TRUE)
   }
   else if (fullPost & !is.null(object.sim$thetaSampList)) {
      if(is.null(thetaSamp)) thetaSamp <- length(object.sim$thetaSampList)
      parRow <- sample(1:nrow(object.sim$thetaSampList[[thetaSamp]]), 1, prob=object.sim$thetaSampList[[thetaSamp]][,1])
      par <- as.vector(object.sim$thetaSampList[[thetaSamp]][parRow,-c(1:3)])
      #print(parRow)
   }
   else par <- object.sim$par
   
   ###
   ### Process parameters for Fortran
   ###
   if (!is.null(err.mfX)) {
      theta.errX <- par[1:n.errX]
      tau2x <- exp(2 * err.mfX %*% theta.errX)
   } else tau2x <- rep(0.0, N)
   if (!is.null(err.mfY)) {
      theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
      tau2y <- exp(2 * err.mfY %*% theta.errY)
   } else tau2y <- tau2x
   theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
   sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
   b <- exp(mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)])
   stay <- rep(0, N)
   if (stopMod) {
      stop.mf <- stop.mf
      theta.stop <- par[(n.errX + n.errY + 2 * n.mov + 1)]
      b <- b / ((1 - stop.mf) ^ exp(theta.stop))
      stay <- ifelse(b==Inf, 1, 0)
      b <- ifelse(b==Inf, 9999, b) 
   }
   if (driftMod) {
      theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                                    (n.errX + n.errY + 2 * n.mov + 2)]
      b.drift <- exp(log(b) - log(1+exp(theta.drift[2])))
      sig2.drift <- exp(log(sig2) + 2 * theta.drift[1])
      call.simulate <- "crwdrift_simulate"
   } else {
      b.drift <- sig2.drift <- 0.0
      call.simulate <- "crw_simulate"
   }
   loctype <- object.sim$loctype
   delta <- as.double(object.sim$delta)
   movMats <- getQT(sig2, b, sig2.drift, b.drift, delta, driftMod)
#   print(movMats$Qmat[1:10,])
   out <- .Fortran(call.simulate,
                   tau2y     = as.double(tau2y),
                   tau2x     = as.double(tau2x),
                   Qmat      = as.double(as.matrix(movMats$Qmat)),
                   Tmat      = as.double(as.matrix(movMats$Tmat)),
                   x         = as.double(object.sim$x),
                   y         = as.double(object.sim$y),
                   loctype   = as.integer(loctype),
                   stay      = as.integer(stay),
                   a1y       = as.double(object.sim$a1.y),
                   a1x       = as.double(object.sim$a1.x),
                   P1y       = as.double(object.sim$P1.y),
                   P1x       = as.double(object.sim$P1.x),
                   lonadj    = as.double(object.sim$lonAdj),
                   N         = as.integer(N),
                   lly       = as.double(0),
                   llx       = as.double(0),
                   alphaY    = as.double(array(rnorm((2+driftMod)*(N+1)),c(2+driftMod,1,N+1))),
                   alphaX    = as.double(array(rnorm((2+driftMod)*(N+1)),c(2+driftMod,1,N+1))),
                   ysim      = as.double(ifelse(loctype==1,9999,rnorm(sum(1-loctype)))),
                   xsim      = as.double(ifelse(loctype==1,9999,rnorm(sum(1-loctype)))),
                   alphaSimY = as.double(matrix(0,N,2+driftMod)),
                   alphaSimX = as.double(matrix(0,N,2+driftMod)),
                   package='crawl')
   if(driftMod) nms <- c("mu","theta","gamma")
   else nms <- c("mu","nu")
   alpha.sim.y <- matrix(out$alphaSimY,N,2+driftMod)
   colnames(alpha.sim.y) <- nms
   alpha.sim.x <- matrix(out$alphaSimX,N,2+driftMod)
   colnames(alpha.sim.x) <- nms
   isw <- ifelse(is.null(object.sim$thetaSampList) & fullPost==TRUE, out$lly+out$llx - object.sim$loglik - dens, 0)
   samp <- list(alpha.sim.y=alpha.sim.y, alpha.sim.x=alpha.sim.x,
                predType=object.sim$predType, Time=object.sim$Time,
                loglik=out$lly+out$llx, par=par, log.isw = isw)
   class(samp) <- c("crwIS","list")
   attr(samp,"coord") <- c(x=object.sim$coord[1], y=object.sim$coord[2])
   attr(samp,"random.drift") <- object.sim$driftMod
   attr(samp,"dry.stop") <- object.sim$stopMod
   attr(samp,"polar.coord") <- object.sim$polar.coord
    if(all(is.nan(alpha.sim.y) | is.nan(alpha.sim.x))) 
    	warning("Simulation failed\n")
   return(samp)
}

