#' Predict animal locations and velocities using a fitted CTCRW model and
#' calculate measurement error fit statistics
#' 

#' 
#' The \code{crwMEfilter} function uses a fitted model object from
#' \code{crwMLE} to predict animal locations (with estimated uncertainty) at
#' times in the original data set and supplimented by times in \code{predTime}.
#' If \code{speedEst} is set to \code{TRUE}, then animal log-speed is also
#' estimated. In addition, the measurement error shock detection filter of de
#' Jong and Penzer (1998) is also calculated to provide a measure for outlier
#' detection.
#' 

#' 
#' The requirements for \code{data} are the same as those for fitting the model
#' in \code{\link{crwMLE}}.
#' 
#' @param object.crwFit A model object from \code{\link{crwMLE}}.
#' @param predTime vector of additional prediction times (numeric or POSIXct). Alternatively, a character vector specifying a time interval (see Details).
#' @param flat logical. Should the result be returned as a flat data.frame.
#' @param ... Additional arguments for testing new features
#' 
#' @details 
#' \itemize{
#' \item("predTime"){
#' \code{predTime} can be either passed as a separate vector of POSIXct or numeric values for additional prediction times beyond the observed times. If the original data were provided as a POSIXct type, then \code{crwPredict} can derive a sequence of regularly spaced prediction times from the original data. This is specified by providing a character string that corresponds to the \code{by} argument of the \code{seq.POSIXt} function (e.g. '1 hour', '30 mins'). \code{crwPredict} will round the first observed time up to the nearest unit (e.g. '1 hour' will round up to the nearest hour, '30 mins' will round up to the nearest minute) and start the sequence from there. The last observation time is truncated down to the nearest unit to specify the end time.
#' }
#' }
#' 
#' @return
#' 
#' List with the following elements:
#' 
#' \item{originalData}{A data.frame with is \code{data} merged with
#' \code{predTime}.}
#' 
#' \item{alpha.hat}{Predicted state}
#' 
#' \item{Var.hat}{array where \code{Var.hat[,,i]} is the prediction
#' covariance matrix for \code{alpha.hat[,i]}.}
#' 
#' \item{fit.test}{A data.frame of chi-square fit (df=2) statistics and naive
#' (pointwise) p-values.}
#' 
#' If \code{flat} is set to \code{TRUE} then a data set is returned with the
#' columns of the original data plus the state estimates, standard errors (se),
#' speed estimates, and the fit statistics and naive p-values.
#' 
#' 
#' @author Devin S. Johnson
#' @references de Jong, P. and Penzer, J. (1998) Diagnosing shocks in time
#' series. Journal of the American Statistical Association 93:796-806.
#' @export

crwPredict=function(object.crwFit, predTime=NULL, flat=TRUE, ...)
{
  if(!exists("getUseAvail")) getUseAvail=FALSE
  if(flat & getUseAvail){
    warning("The 'flat=TRUE' argument cannot be used in conjunction with 'getUseAvail=TRUE' argument.")
    flat <- FALSE
  }
  
  if(inherits(predTime,"character")) {
    t_int <- unlist(strsplit(predTime, " "))
    if(t_int[2] %in% c("min","mins","hour","hours","day","days")) {
      min_dt <- crawl::intToPOSIX(min(object.crwFit$data$TimeNum,na.rm=TRUE))
      max_dt <- crawl::intToPOSIX(max(object.crwFit$data$TimeNum,na.rm=TRUE))
      min_dt <- round(min_dt,t_int[2])
      max_dt <- trunc(max_dt,t_int[2])
      predTime <- seq(min_dt, max_dt, by = predTime)
    } else {
      stop("predTime not specified correctly. see documentation for seq.POSIXt")
    }
  }
  
  ## Model definition/parameters ##
  data <- object.crwFit$data
  driftMod <- object.crwFit$random.drift
  mov.mf <- object.crwFit$mov.mf
  activity <- object.crwFit$activity
  err.mfX <- object.crwFit$err.mfX
  err.mfY <- object.crwFit$err.mfY
  rho = object.crwFit$rho
  par <- object.crwFit$par
  n.errX <- object.crwFit$n.errX
  n.errY <- object.crwFit$n.errY
  n.mov <- object.crwFit$n.mov
  tn <- object.crwFit$Time.name
  if(inherits(predTime, "POSIXct")) predTime <- as.numeric(predTime)#/3600
  
  ## Data setup ##
  if (!is.null(predTime)) {
    if(min(predTime) <  data[1, tn]) {
      warning("Predictions times given before first observation!\nOnly those after first observation will be used.")
      predTime <- predTime[predTime>=data[1,tn]]
    }
    origTime <- data[, tn]
    if (is.null(data$locType)) data$locType <- "o"
    predData <- data.frame(predTime, "p")
    names(predData) <- c(tn, "locType")
    data <- merge(data, predData,
                  by=c(tn, "locType"), all=TRUE)
    dups <- duplicated(data[, tn]) #& data[,"locType"]==1
    data <- data[!dups, ]
    mov.mf <- as.matrix(expandPred(x=mov.mf, Time=origTime, predTime=predTime))
    if (!is.null(activity)) activity <- as.matrix(expandPred(x=activity, Time=origTime, predTime=predTime))
    if (!is.null(err.mfX)) err.mfX <- as.matrix(expandPred(x=err.mfX, Time=origTime, predTime=predTime))
    if (!is.null(err.mfY)) err.mfY <- as.matrix(expandPred(x=err.mfY, Time=origTime, predTime=predTime))
    if (!is.null(rho)) rho <- as.matrix(expandPred(x=rho, Time=origTime, predTime=predTime))
  }
  data$locType[data[,tn]%in%predTime] <- 'p'
  delta <- c(diff(data[, tn]), 1)
  a = object.crwFit$initial.state$a
  P = object.crwFit$initial.state$P
  y = as.matrix(data[,object.crwFit$coord])
  noObs <- as.numeric(is.na(y[,1]) | is.na(y[,2]))
  y[noObs==1,] = 0
  N = nrow(y)
  
  
  
  ###
  ### Process parameters for C++
  ###
  if (!is.null(err.mfX)) {
    theta.errX <- par[1:n.errX]
    Hmat <- exp(2 * err.mfX %*% theta.errX)
  } else Hmat <- rep(0.0, N)
  if (!is.null(err.mfY)) {
    theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
    Hmat <- cbind(Hmat,exp(2 * err.mfY %*% theta.errY))
  } else Hmat <- cbind(Hmat, Hmat)
  if(!is.null(rho)){
    Hmat = cbind(Hmat, sqrt(Hmat[,1])*sqrt(Hmat[,2])*rho)
  } else {Hmat = cbind(Hmat, rep(0,N))}
  Hmat[noObs==1,] = 0
  theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
  sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
  b <- exp(mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)])
  if (!is.null(activity)) {
    theta.activ <- par[(n.errX + n.errY + 2 * n.mov + 1)]
    b <- b / ((activity) ^ exp(theta.activ))
    active <- ifelse(b==Inf, 0, 1)
    b <- ifelse(b==Inf, 0, b) 
  } else active = rep(1,N)
  if (driftMod) {
    theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                         (n.errX + n.errY + 2 * n.mov + 2)]
    b.drift <- exp(log(b) - log(1+exp(theta.drift[2])))
    sig2.drift <- exp(log(sig2) + 2 * theta.drift[1])
    out = CTCRWPREDICT_DRIFT(y, Hmat, b, b.drift, sig2, sig2.drift, delta, noObs, active, a, P)
  } else {
    out=CTCRWPREDICT(y, Hmat, b, sig2, delta, noObs, active, a, P)
  }
  
  pred <- data.frame(t(out$pred))
  if (driftMod) {
    names(pred) <- c("mu.x", "theta.x", "gamma.x","mu.y", "theta.y", "gamma.y")
  } else names(pred) <- c("mu.x", "nu.x", "mu.y","nu.y")
  var <- zapsmall(out$predVar)
  obsFit <- data.frame(predObs.x=out$predObs[1,],
                       predObs.y=out$predObs[2,])
  obsFit$outlier.chisq <- as.vector(out$chisq)
  obsFit$naive.p.val <- 1 - pchisq(obsFit$outlier.chisq, 2)
  if(getUseAvail){
    warning("'getUseAvail' not implemented yet in this version of 'crawl' contact maintainer to fix this! ")
    #     idx <- data$locType=="p"
    #     movMatsPred <- getQT(sig2[idx], b[idx], sig2.drift[idx], b.drift[idx], delta=c(diff(data[idx,tn]),1), driftMod)
    #     TmatP <- movMatsPred$Tmat
    #     QmatP <- movMatsPred$Qmat
    #     avail <- t(sapply(1:(nrow(TmatP)-1), makeAvail, Tmat=TmatP, Qmat=QmatP, predx=predx[idx,], predy=predy[idx,], 
    #                       vary=vary[,,idx], varx=varx[,,idx], driftMod=driftMod, lonadj=lonAdjVals[idx]))
    #     avail <- cbind(data[idx,tn][-1], avail)
    #     colnames(avail) <- c(tn, "meanAvail.x", "meanAvail.y", "varAvail.x", "varAvail.y")
    #     use <- cbind(data[idx,tn], predx[idx,1], predy[idx,1], varx[1,1,idx], vary[1,1,idx])[-1,]
    #     colnames(use) <- c(tn, "meanUse.x", "meanUse.y", "varUse.x", "varUse.y")
    #     UseAvail.lst <- list(use=use, avail=avail)
    UseAvail.lst=NULL
  }
  else UseAvail.lst=NULL
  speed = sqrt(apply(as.matrix(pred[,2:(2+driftMod)]), 1, sum)^2 + 
                 apply(as.matrix(pred[,(4+driftMod):(4+2*driftMod)]), 1, sum)^2)
  out <- list(originalData=fillCols(data), alpha.hat=pred, 
              V.hat=var, speed=speed, loglik=out$ll, useAvail=UseAvail.lst)
  if (flat) {
    out <- cbind(fillCols(crawl::flatten(out)), obsFit)
    attr(out, "flat") <- TRUE
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "activity.model") <- !is.null(object.crwFit$activity)
    attr(out, "Time.name") <- tn
  } else {
    out <- append(out, list(fit.test=obsFit))
    attr(out, "flat") <- FALSE
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "activity.model") <- !is.null(object.crwFit$activity)
    attr(out, "Time.name") <- tn
  }
  class(out) <- c(class(out),"crwPredict")
  return(out)
}
