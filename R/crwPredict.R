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
  fixPar = object.crwFit$fixPar
  theta = object.crwFit$par[is.na(fixPar)]
  
  
  ## the typical expectation is for tn to be POSIXct. But, some users may decide
  ## to pass a numeric time vector. Here, we'll confirm numeric or POSIXct and
  ## set return_posix to TRUE unless the submitted values are numeric. In the
  ## case where the data values are one and the predTime values are another, we
  ## will convert to POSIXct and return POSIXct
  ## 
  return_posix <- ifelse((inherits(predTime,"POSIXct") | inherits(predTime, "character")) & 
                           inherits(data[,tn],"POSIXct"), 
                         TRUE, FALSE)
  if(!return_posix) {
    # if(inherits(predTime,"numeric") && inherits(data[, tn],"numeric")) {
    #   warning("numeric time values detected. numeric values will be returned.")
    # } 
    if(inherits(predTime,"numeric") && inherits(data[, tn], "POSIXct")) {
      # warning("predTime provided as numeric. converting it to POSIXct.")
      stop("predTime provided as numeric and original time data was POSIX!")
      # predTime <- lubridate::as_datetime(predTime)
    }
    if(inherits(predTime,"POSIXct") && inherits(data[, tn], "numeric")) {
      # warning("input data time column provided as numeric. converting to POSIXct")
      stop("predTime provided as POSIX and original data was numeric!")
      # data[, tn] <- lubridate::as_datetime(data[, tn])
    }
  }
  
  if(!is.null(predTime)){
    
    if(inherits(predTime,"character")) {
      if(!inherits(data[,tn],"POSIXct")) stop("Character specification of predTime can only be used with POSIX times in the original data!")
      t_int <- unlist(strsplit(predTime, " "))
      if(t_int[2] %in% c("min","mins","hour","hours","day","days")) {
        min_dt <- min(data[,tn],na.rm=TRUE)
        max_dt <- max(data[,tn],na.rm=TRUE)
        min_dt <- lubridate::ceiling_date(min_dt,t_int[2])
        max_dt <- lubridate::floor_date(max_dt,t_int[2])
        predTime <- seq(min_dt, max_dt, by = predTime)
      } else {
        stop("predTime not specified correctly. see documentation for seq.POSIXt")
      }
    }
    
    if(inherits(predTime, "POSIXct")){
      ts = attr(object.crwFit, "time.scale")
      predTime = as.numeric(predTime)/ts
    }
    
    ## Data setup ##
    if(min(predTime) <  min(data$TimeNum)) {
      warning("Predictions times given before first observation!\nOnly those after first observation will be used.")
      predTime <- predTime[predTime>=data$TimeNum]
    }
    origTime <- data$TimeNum
    if (is.null(data$locType)) {
      data$locType <- "o"
    }
    predData <- data.frame(predTime, "p")
    names(predData) <- c("TimeNum", "locType")
    # predTime <- as.numeric(predTime)
    data <- merge(data, predData,
                  by=c("TimeNum", "locType"), all=TRUE)
    dups <- duplicated(data$TimeNum) #& data[,"locType"]==1
    data <- data[!dups, ]
    mov.mf <- as.matrix(expandPred(x=mov.mf, Time=origTime, predTime=predTime))
    if (!is.null(activity)) activity <- as.matrix(expandPred(x=activity, Time=origTime, predTime=predTime))
    if (!is.null(err.mfX)) err.mfX <- as.matrix(expandPred(x=err.mfX, Time=origTime, predTime=predTime))
    if (!is.null(err.mfY)) err.mfY <- as.matrix(expandPred(x=err.mfY, Time=origTime, predTime=predTime))
    if (!is.null(rho)) rho <- as.matrix(expandPred(x=rho, Time=origTime, predTime=predTime))
  }
  data$locType[data$TimeNum%in%predTime] <- 'p'
  delta <- c(diff(data$TimeNum), 1)
  y = as.matrix(data[,object.crwFit$coord])
  noObs <- as.numeric(is.na(y[,1]) | is.na(y[,2]))
  y[noObs==1,] = 0
  N = nrow(y)
  
  ###
  ### Process parameters for C++
  ###
  argslist = par2arglist(theta, fixPar, y, noObs, delta,
                         mov.mf, err.mfX, err.mfY, rho, activity,
                         n.errX, n.errY, n.mov, driftMod)
  if (driftMod) {
    out = CTCRWPREDICT_DRIFT(y, argslist$Hmat, argslist$b, argslist$b.drift, argslist$sig2, 
                             argslist$sig2.drift, delta, noObs, argslist$active, argslist$a, argslist$P)
  } else {
    out=CTCRWPREDICT(y, argslist$Hmat, argslist$b, argslist$sig2, delta, noObs, argslist$active, argslist$a, argslist$P)
  }
  
  pred <- data.frame(t(out$pred))
  if (driftMod) {
    names(pred) <- c("mu.x", "theta.x", "gamma.x","mu.y", "theta.y", "gamma.y")
  } else names(pred) <- c("mu.x", "nu.x", "mu.y","nu.y")
  var <- zapsmall(out$predVar)
  
  speed = sqrt(apply(as.matrix(pred[,2:(2+driftMod)]), 1, sum)^2 + 
                 apply(as.matrix(pred[,(4+driftMod):(4+2*driftMod)]), 1, sum)^2)
  
  obsFit <- data.frame(predObs.x=out$predObs[1,],
                       predObs.y=out$predObs[2,])
  obsFit$outlier.chisq <- as.vector(out$chisq)
  obsFit$naive.p.val <- 1 - pchisq(obsFit$outlier.chisq, 2)
  
  out <- list(originalData=fillCols(data), alpha.hat=pred, 
              V.hat=var, speed=speed, loglik=out$ll)
  if(return_posix){
    out$originalData[,tn] = lubridate::as_datetime(out$originalData$TimeNum*ts)
  } else out$originalData[,tn] = out$originalData$TimeNum
  
  # if(getUseAvail){
  #   idx <- data$locType=="p"
  #   movMatsPred <- getQT(sig2[idx], b[idx], sig2.drift[idx], b.drift[idx], delta=c(diff(data[idx,tn]),1), driftMod)
  #   TmatP <- movMatsPred$Tmat
  #   QmatP <- movMatsPred$Qmat
  #   avail <- t(sapply(1:(nrow(TmatP)-1), makeAvail, Tmat=TmatP, Qmat=QmatP, predx=predx[idx,], predy=predy[idx,], 
  #                     vary=vary[,,idx], varx=varx[,,idx], driftMod=driftMod, lonadj=lonAdjVals[idx]))
  #   avail <- cbind(data[idx,tn][-1], avail)
  #   colnames(avail) <- c(tn, "meanAvail.x", "meanAvail.y", "varAvail.x", "varAvail.y")
  #   use <- cbind(data[idx,tn], predx[idx,1], predy[idx,1], varx[1,1,idx], vary[1,1,idx])[-1,]
  #   colnames(use) <- c(tn, "meanUse.x", "meanUse.y", "varUse.x", "varUse.y")
  #   UseAvail.lst <- list(use=use, avail=avail)
  # }
  # else UseAvail.lst=NULL
  
  if (return.type == "flat") {
    out <- fillCols(crawl::flatten(out))
    attr(out, "flat") <- TRUE
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "activity.model") <- !is.null(object.crwFit$activity)
    attr(out, "Time.name") <- tn
    attr(out, "time.scale") = ts
    attr(out,"epsg") <- attr(object.crwFit,"epsg")
    attr(out,"proj4") <- attr(object.crwFit,"proj4")
  } else if (return.type == "list") {
    out <- append(out, list(fit.test=obsFit))
    attr(out, "flat") <- FALSE
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "activity.model") <- !is.null(object.crwFit$activity)
    attr(out, "Time.name") <- tn
    attr(out, "time.scale") = ts
    attr(out,"epsg") <- attr(object.crwFit,"epsg")
    attr(out,"proj4") <- attr(object.crwFit,"proj4")
  } else if (return.type == "minimal") {
    out <- fillCols(out$originalData)
    out <- cbind(out, pred)
    attr(out, "flat") <- TRUE
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "activity.model") <- !is.null(object.crwFit$activity)
    attr(out, "Time.name") <- tn
    attr(out, "time.scale") = ts
    attr(out,"epsg") <- attr(object.crwFit,"epsg")
    attr(out,"proj4") <- attr(object.crwFit,"proj4")
  }
  class(out) <- c(class(out),"crwPredict")
  return(out)
}
