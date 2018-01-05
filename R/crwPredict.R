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
#' @param predTime vector of desired prediction times (numeric or POSIXct). Alternatively, a character vector specifying a time interval (see Details).
#' @param return.type character. Should be one of \code{"minimal","flat","list"} (see Details).
#' @param ... Additional arguments for testing new features
#' 
#' @details 
#' \itemize{
#' \item("predTime"){
#' \code{predTime} can be either passed as a separate vector of POSIXct or numeric values for all prediction times expected in the returned object. Note, previous versions of \code{crwPredict} would return both times specified via \code{predTime} as well as each original observed time. This is no longer the default (see \item{return.type}). If the original data were provided as a POSIXct type, then \code{crwPredict} can derive a sequence of regularly spaced prediction times from the original data. This is specified by providing a character string that corresponds to the \code{by} argument of the \code{seq.POSIXt} function (e.g. '1 hour', '30 mins'). \code{crwPredict} will round the first observed time up to the nearest unit (e.g. '1 hour' will round up to the nearest hour, '30 mins' will round up to the nearest minute) and start the sequence from there. The last observation time is truncated down to the nearest unit to specify the end time.
#' }
#' }
#' 
#' @return
#' 
#' There are three possible return types specified with \code{return.type}:
#' 
#' \item{minimal}{a data.frame with a minimal set of columns: 
#' \code{date_time,mu.x,mu.y,se.mu.x,se.mu.y}}
#' 
#' \item{flat}{a data set is returned with the
#' columns of the original data plus the state estimates, standard errors (se),
#' and speed estimates}
#' 
#' \item{list}{List with the following elements:}
#' 
#' \item{originalData}{A data.frame with \code{data} merged with
#' \code{predTime}.}
#' 
#' \item{alpha.hat}{Predicted state}
#' 
#' \item{Var.hat}{array where \code{Var.hat[,,i]} is the prediction
#' covariance matrix for \code{alpha.hat[,i]}.}
#' 
#' 
#' @author Devin S. Johnson
#' @references de Jong, P. and Penzer, J. (1998) Diagnosing shocks in time
#' series. Journal of the American Statistical Association 93:796-806.
#' @export

crwPredict=function(object.crwFit, predTime=NULL, return.type="minimal", ...)
{
  data <- object.crwFit$data
  tn <- object.crwFit$Time.name
  
  ## the typical expectation is for tn to be POSIXct. But, some users may decide
  ## to pass a numeric time vector. Here, we'll confirm numeric or POSIXct and
  ## set return_posix to TRUE unless the submitted values are numeric. In the
  ## case where the data values are one and the predTime values are another, we
  ## will convert to POSIXct and return POSIXct
  ## 
  return_posix <- ifelse(inherits(predTime,"POSIXct") && 
                           inherits(data[,tn],"POSIXct"), 
                         TRUE, FALSE)
  if(!return_posix) {
    if(inherits(predTime,"numeric") && inherits(data[, tn],"numeric")) {
      message("numeric time values detected. numeric values will be returned.")
    } 
    if(inherits(predTime,"numeric") && inherits(data[, tn], "POSIXct")) {
      message("predTime provided as numeric. converting to POSIXct.")
      predTime <- lubridate::as_datetime(predTime)
    }
    if(inherits(predTime,"POSIXct") && inherits(data[, data], "numeric")) {
      message("input data time column provided as numeric. converting to POSIXct")
      data[, tn] <- lubridate::as_datetime(data[, tn])
    }
  }

  
  if(inherits(predTime,"character")) {
    if(inherits(data[, tn], "numeric")) {
      warning("predTime specified as character string and data time column as numeric. converting data time column to POSIXct.")
      data[, tn] <- lubridate::as_datetime(data[, tn])
    }
    t_int <- unlist(strsplit(predTime, " "))
    if(t_int[2] %in% c("min","mins","hour","hours","day","days")) {
      if(!inherits(data[tn],"POSIXct")) {
      min_dt <- min(data[tn],na.rm=TRUE)
      max_dt <- max(data[tn],na.rm=TRUE)
      }
      min_dt <- round(min_dt,t_int[2])
      max_dt <- trunc(max_dt,t_int[2])
      predTime <- seq(min_dt, max_dt, by = predTime)
    } else {
      stop("predTime not specified correctly. see documentation for seq.POSIXt")
    }
  }
  
  ## Model definition/parameters ##
  
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
  
  ## Data setup ##
  if (!is.null(predTime)) {
    if(min(predTime) <  data[1, tn]) {
      warning("Predictions times given before first observation!\nOnly those after first observation will be used.")
      predTime <- predTime[predTime>=data[1,tn]]
    }
    origTime <- as.numeric(data[, tn])
    predTime <- as.numeric(predTime)
    if (is.null(data$locType)) {
      data$locType <- "o"
    }
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
  
  speed = sqrt(apply(as.matrix(pred[,2:(2+driftMod)]), 1, sum)^2 + 
                 apply(as.matrix(pred[,(4+driftMod):(4+2*driftMod)]), 1, sum)^2)
  out <- list(originalData=fillCols(data), alpha.hat=pred, 
              V.hat=var, speed=speed, loglik=out$ll)
  if (return.type == "flat") {
    out <- fillCols(crawl::flatten(out))
    attr(out, "flat") <- TRUE
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "activity.model") <- !is.null(object.crwFit$activity)
    attr(out, "Time.name") <- tn
  } else if (return.type == "list") {
    attr(out, "flat") <- FALSE
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "activity.model") <- !is.null(object.crwFit$activity)
    attr(out, "Time.name") <- tn
  } else if (return.type == "minimal") {
    out <- fillCols(data)[,c(tn,"locType")]
    #out <- cbind(var[1,1,], var[1,3,])
    
    #out <- out[,c(tn,"mu.x","mu.y","se.mu.x","se.mu.y")]
  }
  if (return_posix) {
    out$tn <- lubridate::as_datetime(out$tn)
  }
  class(out) <- c(class(out),"crwPredict")
  return(out)
}
