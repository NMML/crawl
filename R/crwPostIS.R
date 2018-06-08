#' Simulate a value from the posterior distribution of a CTCRW model
#' 

#' 
#' The crwPostIS draws a set of states from the posterior distribution of a
#' fitted CTCRW model. The draw is either conditioned on the fitted parameter
#' values or "full" posterior draw with approximated parameter posterior
#' 

#' 
#' The crwPostIS draws a posterior sample of the track state matrices. If
#' fullPost was set to TRUE when the object.sim was build in
#' \link{crwSimulator} then a psuedo-posterior draw will be made by first
#' sampling a parameter value from a multivariate t distribution which
#' approximates the marginal posterior distribution of the parameters. The
#' covariance matrix from the fitted model object is used to scale the MVt
#' approximation. In addition, the factor "scale" can be used to further adjust
#' the approximation. Further, the parameter simulations are centered on the
#' fitted values.
#' 
#' To correct for the MVt approximation, the importance sampling weight is also
#' supplied. When calulating averages of track functions for Bayes estimates
#' one should use the importance sampling weights to calculate a weighted
#' average (normalizing first, so the weights sum to 1).
#' 
#' @param object.sim A crwSimulator object from \code{\link{crwSimulator}}.
#' @param fullPost logical. Draw parameter values as well to simulate full
#' posterior
#' @param df degrees of freedom for multivariate t distribution approximation
#' to parameter posterior
#' @param scale Extra scaling factor for t distribution approximation
#' @param thetaSamp If multiple parameter samples are available in object.sim,
#' setting \code{thetaSamp=n} will use the nth sample. Defaults to the last.
#' @return
#' 
#' List with the following elements:
#' 
#' \item{alpha.sim.y}{A matrix a simulated latitude state values}
#' 
#' \item{alpha.sim.x}{Matrix of simulated longitude state values}
#' 
#' \item{locType}{Indicates prediction types with a "p" or observation times
#' with an "o"} \item{Time}{Initial state covariance for latitude}
#' 
#' \item{loglik}{log likelihood of simulated parameter}
#' 
#' \item{par}{Simulated parameter value}
#' 
#' \item{log.isw}{non normalized log importance sampling weight}
#' @author Devin S. Johnson
#' @seealso See \code{demo(northernFurSealDemo)} for example.
#' @export

crwPostIS = function(object.sim, fullPost=TRUE, df=Inf, scale=1, thetaSamp=NULL)
    ################################################################################
################################################################################
{
  if(!inherits(object.sim, 'crwSimulator')) stop("Argument needs to be of class 'crwSimulator'\nUse 'crwSimulator( )' to create")
  fixPar <- object.sim$fixPar
  Cmat <- object.sim$Cmat[is.na(fixPar),is.na(fixPar)]
  se <- sqrt(diag(Cmat))
  err.mfX <- object.sim$err.mfX
  err.mfY <- object.sim$err.mfY
  par <- object.sim$par
  n2ll.mode <- -2*object.sim$loglik
  activity <- object.sim$activity
  driftMod <- object.sim$driftMod
  mov.mf <- object.sim$mov.mf
  y <- object.sim$y
  noObs = object.sim$noObs
  delta <- object.sim$delta
  n.errX <- object.sim$n.errX
  n.errY <- object.sim$n.errY
  rho = object.sim$rho
  n.mov <- object.sim$n.mov
  N <- object.sim$N
  lower <- object.sim$lower
  upper <- object.sim$upper 
  prior <- object.sim$prior
  eInd <- is.na(fixPar)
  ###
  ### Sample parameter vector
  ###
  if(fullPost){
    if(is.null(object.sim$thetaSampList)){
      eps <- rmvtt(mu=rep(0,sum(eInd)), Sigma=scale*Cmat, df=df, lower-par[eInd], upper-par[eInd])
      par[eInd] <- par[eInd] + eps
      if(df==Inf) dens <- dmvnorm(eps, sigma=scale*Cmat, log=TRUE) - dmvnorm(0.0*eps, sigma=scale*Cmat, log=TRUE)
      else dens <- dmvt(eps, sigma=scale*Cmat, df=df, log=TRUE) - dmvt(0.0*eps, sigma=scale*Cmat, df=df, log=TRUE)
    } else{
      if(is.null(thetaSamp)) thetaSamp <- length(object.sim$thetaSampList)
      parRow <- sample(1:nrow(object.sim$thetaSampList[[thetaSamp]]), 1, prob=object.sim$thetaSampList[[thetaSamp]][,1])
      par <- as.vector(object.sim$thetaSampList[[thetaSamp]][parRow,-c(1:3)])
    }
  }
  
  ###
  ### Process parameters for C++
  ###
  theta = object.sim$par[is.na(object.sim$fixPar)]
  argslist = par2arglist(theta, fixPar, y, noObs, delta,
                         mov.mf, err.mfX, err.mfY, rho, activity,
                         n.errX, n.errY, n.mov, driftMod)
  
  if (driftMod) {
    out=CTCRWSAMPLE_DRIFT(y, argslist$Hmat,  argslist$b,  argslist$b.drift,  argslist$sig2, 
                          argslist$sig2.drift, delta, noObs,  argslist$active,  argslist$a,  argslist$P)
  } else {
    out=CTCRWSAMPLE(y,  argslist$Hmat,  argslist$b,  argslist$sig2,  delta, noObs,  argslist$active,  argslist$a,  argslist$P)
  }
  
  if(driftMod){
    colnames(out$sim) <- apply(expand.grid(c("mu","theta","gamma"), c("x","y")), 1, paste, collapse=".")
  }  else {
    colnames(out$sim) <- apply(expand.grid(c("mu","nu"), c("x","y")), 1, paste, collapse=".")
  }
  ln.prior = ifelse(!is.null(object.sim$prior), object.sim$prior(par[eInd]), 0)
  isw <- ifelse(is.null(object.sim$thetaSampList) & fullPost==TRUE, out$ll - object.sim$loglik - dens, 0) + ln.prior
  samp <- list(alpha.sim=out$sim,
               locType=object.sim$locType, TimeNum=object.sim$TimeNum, 
               loglik=out$lly+out$llx, par=par, log.isw = isw)
  samp[[object.sim$Time.name]] = object.sim$TimeNum
  if(object.sim$return_posix) samp[[object.sim$Time.name]] = lubridate::as_datetime(samp[[object.sim$Time.name]])
  class(samp) <- c("crwIS","list")
  attr(samp, "Time.name") = object.sim$Time.name
  attr(samp,"coord") <- object.sim$coord
  attr(samp,"random.drift") <- object.sim$driftMod
  attr(samp,"activity.model") <- !is.null(object.sim$activity)
  attr(samp,"epsg") <- attr(object.sim,"epsg")
  attr(samp,"proj4") <- attr(object.sim,"proj4")
  return(samp)
}

