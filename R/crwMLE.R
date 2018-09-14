#' Fit Continuous-Time Correlated Random Walk Models to Animal Telemetry Data
#' 
#' The function uses the Kalman filter to estimate movement paramters in a
#' state-space version of the continuous-time movement model. Separate models
#' are specified for movement portion and the location error portion. Each
#' model can depend on time indexed covariates. A \dQuote{haul out} model where
#' movement is allowed to completely stop, as well as, a random drift model can
#' be fit with this function.
#' 

#' 
#' A full model specification involves 4 components: a movement model, an
#' activity model, 2 location error models, and a drift indication. The
#' movement model (\code{mov.model}) specifies how the movement parameters
#' should vary over time. This is a function of specified, time-indexed,
#' covariates. The movement parameters (sigma for velocity variation and beta
#' for velocity autocorrelation) are both modeled with a log link as par =
#' exp(eta), where eta is the linear predictor based on the covariates. The
#' \code{err.model} specification is a list of 2 such models, one for
#' \dQuote{longitude} and one for \dQuote{latitude} (in that order) location
#' error. If only one location error model is given, it is used for both
#' coordinates (parameter values as well). If \code{drift.model} is set to
#' \code{TRUE}, then, 2 additional parameters are estimated for the drift
#' process, a drift variance and a beta multiplier. 
#' 
#' \code{theta} and \code{fixPar} are vectors with the appropriate number or
#' parameters. \code{theta} contains only those paraemters which are to be
#' estimated, while \code{fixPar} contains all parameter values with \code{NA}
#' for parameters which are to be estimated.
#' 
#' The data set specified by \code{data} must contain a numeric or POSIXct column which is
#' used as the time index for analysis. The column name is specified by the
#' \code{Time.name} argument. If a POSIXct column is used it is internally converted to a
#' numeric vector with units of \code{time.scale}. Also, for activity models, the
#' sactivity covariate must be between 0 and 1 inclusive, with 0 representing complete stop
#' of the animal (no true movement, however, location error can still occur) and 1 
#' represent unhindered movement. The coordinate location should have \code{NA} where no
#' location is recorded, but there is a change in the movment covariates.
#' 
#' The CTCRW models can be difficult to provide good initial values for
#' optimization. If \code{initialSANN} is specified then simulated annealing is
#' used first to obtain starting values for the specified optimaization method.
#' If simulated annealing is used first, then the returned \code{init} list of
#' the crwFit object will be a list with the results of the simulated annealing
#' optimization.
#' 
#' @param mov.model formula object specifying the time indexed covariates for
#' movement parameters.
#' @param err.model A 2-element list of formula objects specifying the time
#' indexed covariates for location error parameters.
#' @param activity formula object giving the covariate for the activity (i.e., stopped or fully moving)
#' portion of the model.
#' @param drift logical indicating whether or not to include a random
#' drift component. For most data this is usually not necessary. See \code{\link{northernFurSeal}} for an example
#' using a drift model.
#' @param data data.frame object containg telemetry and covariate data. A 
#'   'SpatialPointsDataFrame' object from the package 'sp' or an 'sf' object
#'   from the 'sf' package with a geometry column of type \code{sfc_POINT}.
#'   'spacetime' objects were previously accepted but no longer valid. 
#'   Values for coords will be taken from 
#'   the spatial data set and ignored in the arguments. Spatial data must have a
#'   valid proj4string or epsg and must NOT be in longlat.
#' @param coord A 2-vector of character values giving the names of the "X" and
#' "Y" coordinates in \code{data}.
#' @param Time.name character indicating name of the location time column
#' @param time.scale character. Scale for conversion of POSIX time to numeric for modeling. Defaults to "hours".
#' @param theta starting values for parameter optimization.
#' @param fixPar Values of parameters which are held fixed to the given value.
#' @param method Optimization method that is passed to \code{\link{optim}}.
#' @param control Control list which is passed to \code{\link{optim}}.
#' @param constr Named list with elements \code{lower} and \code{upper} that
#' are vectors the same length as theta giving the box constraints for the
#' parameters
#' @param prior A function returning the log-density function of the parameter
#' prior distribution. THIS MUST BE A FUNCTION OF ONLY THE FREE PARAMETERS. Any 
#' fixed parameters should not be included.
#' @param need.hess A logical value which decides whether or not to evaluate
#' the Hessian for parameter standard errors
#' @param initialSANN Control list for \code{\link{optim}} when simulated
#' annealing is used for obtaining start values. See details
#' @param attempts The number of times likelihood optimization will be
#' attempted
#' @param ... Additional arguments that are ignored.
#' @return
#' 
#' A list with the following elements:
#' 
#' \item{par}{Parameter maximum likelihood estimates (including fixed
#' parameters)}
#' 
#' \item{estPar}{MLE without fixed parameters}
#' 
#' \item{se}{Standard error of MLE}
#' 
#' \item{ci}{95\% confidance intervals for parameters}
#' 
#' \item{Cmat}{Parameter covariance matrix}
#' 
#' \item{loglik}{Maximized log-likelihood value}
#' 
#' \item{aic}{Model AIC value}
#' 
#' \item{coord}{Coordinate names provided for fitting}
#' 
#' \item{fixPar}{Fixed parameter values provided}
#' 
#' \item{convergence}{Indicator of convergence (0 = converged)}
#' 
#' \item{message}{Meesages given by \code{optim} during parameter optimization}
#' 
#' \item{activity}{Model provided for stopping variable}
#' 
#' \item{drift}{Logical value indicating random drift model}
#' 
#' \item{mov.model}{Model description for movement component}
#' 
#' \item{err.model}{Model description for location error component}
#' 
#' \item{n.par}{number of parameters}
#' 
#' \item{nms}{parameter names}
#' 
#' \item{n.mov}{number of movement parameters}
#' 
#' \item{n.errX}{number or location error parameters for "longitude" error
#' model}
#' 
#' \item{n.errY}{number or location error parameters for "latitude" error
#' model}
#' 
#' \item{stop.mf}{covariate for stop indication in stopping models}
#' 
#' \item{polar.coord}{Logical indicating coordinates are polar latitude and
#' longitude}
#' 
#' \item{init}{Initial values for parameter optimization}
#' 
#' \item{data}{Original data.frame used to fit the model}
#' 
#' \item{lower}{The lower parameter bounds}
#' 
#' \item{upper}{The upper parameter bounds}
#' 
#' \item{need.hess}{Logical value}
#' 
#' \item{runTime}{Time used to fit model}
#' @author Devin S. Johnson, Josh M. London
#' @export

crwMLE = function(mov.model=~1, err.model=NULL, activity=NULL, drift=FALSE,
                  data, coord=c("x", "y"), Time.name="time", time.scale="hours", #initial.state, 
                  theta, fixPar, method="Nelder-Mead", control=NULL, constr=list(lower=-Inf, upper=Inf), 
                  prior=NULL, need.hess=TRUE, initialSANN=list(maxit=200), attempts=1, ...)
{
  #if(drift) stop("At this time drift models are not supported with this function. Use 'crwMLE' for now.\n")
  st <- Sys.time()
  
  ### Transform 'sp' package SpatialPointsDataFrame
  if(inherits(data, "trip")){
    Time.name <- data@TOR.columns[1]
  }
  p4 <- list(epsg = NULL, proj4string = NULL)
  if(inherits(data, "SpatialPoints")) {	
    if(!sp::is.projected(data)) {
      stop("proj4string for data of sp class is not specified.")
    }
    if("+proj=longlat" %in% strsplit(sp::proj4string(data), " ")[[1]]) stop("Location data is provided in longlat; must be projected.")	
    p4$proj4string <- sp::proj4string(data)
    coordVals <- as.data.frame(sp::coordinates(data))	
    coord <- names(coordVals)	
    data <- cbind(slot(data,"data"), coordVals)    
  }
  if(inherits(data,"sf") && inherits(sf::st_geometry(data),"sfc_POINT")) {
    if (sf::st_is_longlat(data)) {
      stop("Location data is provided in longlat; must be projected.")
    }
    p4$epsg <- sf::st_crs(data)$epsg
    p4$proj4string <- sf::st_crs(data)$proj4string
    if(!any(names(data) %in% c("x","y"))) {
      # warning("no 'x' and 'y' columns detected in 'sf' object so will create")
      coordVals <- as.data.frame(do.call(rbind,sf::st_geometry(data)))
      coordVals <- stats::setNames(coordVals, c("x","y"))
      sf::st_geometry(data) <- NULL
      data <- cbind(data, coordVals)
    } else {
      warning("'x' and 'y' columns detected in 'sf' so will use these")
      sf::st_geometry(data) <- NULL
    }
  }
  if(inherits(data,"tbl_df")) {
    data <- as.data.frame(data)
  }
  if(inherits(data[,Time.name],"POSIXct")){
    if(time.scale %in% c("hours", "hour")){
      ts = 60*60
    } else if(time.scale %in% c("days", "day")){
      ts = 60*60*24
    } else if(time.scale %in% c("sec","secs","second","seconds")){
      ts = 1
    } else if(time.scale %in% c("min","mins","minute","minutes")){
      ts = 60
    } else stop("'time.scale' not specified correctly!")
    data$TimeNum <- as.numeric(data[,Time.name])/ts
  } else{
    data$TimeNum <- as.numeric(data[,Time.name])
    ts = 1
  }
  
  ## SET UP MODEL MATRICES AND PARAMETERS ##
  errMod <- !is.null(err.model)
  #if(!errMod) stop("Error model must be specified! (argument 'err.model' is currently set to NULL)")
  activeMod <- !is.null(activity)
  driftMod <- drift
  # if (
  #   length(initial.state$a) != 2*(driftMod+2) | all(dim(initial.state$P) != c(2*(driftMod+2), 2*(driftMod+2))) 
  # ) stop("Dimentions of 'initial.state' argument are not correct for the specified model")
  mov.mf <- model.matrix(mov.model, model.frame(mov.model, data, na.action=na.pass))
  if (any(is.na(mov.mf))) stop("Missing values are not allowed in movement covariates!")
  n.mov <- ncol(mov.mf)
  if (errMod) {
    err.mfX <- model.matrix(err.model$x,model.frame(err.model$x, data, na.action=na.pass))
    err.mfX <- ifelse(is.na(err.mfX), 0, err.mfX)
    n.errX <- ncol(err.mfX)
    if (!is.null(err.model$y)) {
      err.mfY <- model.matrix(err.model$y,model.frame(err.model$y, data, na.action=na.pass))
      err.mfY <- ifelse(is.na(err.mfY), 0, err.mfY)
      n.errY <- ncol(err.mfY)
    } else {
      err.mfY <- NULL
      n.errY <- 0
    }
    if(!is.null(err.model$rho)){
      rho = model.matrix(err.model$rho,model.frame(err.model$rho, data, na.action=na.pass))[,-1]
      if(any(rho > 1 | rho < -1, na.rm=TRUE)) stop("Error model correlation outside of the range (-1, 1).")
      rho <- ifelse(is.na(rho), 0, rho)
    } else rho = NULL
  } else {
    n.errY <- n.errX <- 0
    err.mfX <- err.mfY <- rho <- NULL
  }
  if (activeMod) {
    #stop.model
    activity <- model.matrix(activity, model.frame(activity, data, na.action=na.pass))
    if (ncol(activity) > 2) stop("There can only be one activity variable.")
    activity <- as.double(activity[,2])
    if (any(activity < 0) | any(activity > 1)) stop("'activity' variable must be >=0 and <=1.")
    if (any(is.na(activity))) stop("Missing values are not allowed in the activity variable.")
  } else activity <- NULL
  n.drift <- as.integer(driftMod)
  n.activ <- as.integer(activeMod)
  b.nms <- paste("ln beta ", colnames(mov.mf), sep="")
  sig.nms <- paste("ln sigma ", colnames(mov.mf), sep="")
  if (errMod) {
    if (!is.null(err.model$y)) {
      tau.nms <- c(paste("ln tau.x ", colnames(err.mfX), sep=""),
                   paste("ln tau.y ", colnames(err.mfY), sep=""))
    } else tau.nms <- paste("ln tau ", colnames(err.mfX), sep="")
  } else tau.nms <- NULL
  if (activeMod){
    active.nms <- "ln phi"
    } else active.nms <- NULL
  if (driftMod) {
    drift.nms <- c("ln sigma.drift/sigma", "ln psi-1")
  } else drift.nms <- NULL
  nms <- c(tau.nms, sig.nms, b.nms, active.nms, drift.nms)
  n.par <- length(nms)
  if (missing(fixPar)) fixPar <- rep(NA, n.par)
  n.theta = sum(is.na(fixPar))
  if (length(fixPar)!=n.par) stop("'fixPar' argument is not the right length! The number of parameters in the model is ", n.par, "\n")
  if(!(length(constr$lower)==1 | length(constr$lower)==sum(is.na(fixPar)))) stop("The number of lower contraints specified is not correctly! The number of free parameters is ", sum(is.na(fixPar)),"\n")
  if(!(length(constr$upper)==1 | length(constr$upper)==sum(is.na(fixPar)))) stop("The number of upper contraints specified is not correctly! The number of free parameters is ", sum(is.na(fixPar)),"\n")
  # if(length(constr$upper)==1) constr$upper <- rep(constr$upper, sum(is.na(fixPar)))
  # if(length(constr$lower)==1) constr$lower <- rep(constr$lower, sum(is.na(fixPar)))
  if (missing(theta)) theta = rep(0,n.theta)
  theta[theta<constr$lower] = constr$lower[theta<constr$lower] + 0.01
  theta[theta>constr$upper] = constr$upper[theta>constr$upper] - 0.01
  if(driftMod & is.na(fixPar[n.par])) theta[sum(is.na(fixPar))] <- log(diff(range(data$TimeNum))/9)
  if (length(theta) != n.theta) {
    stop("\nWrong number of parameters specified in start value.\n")
  }
  
  y = data[,c(coord[1],coord[2])]
  noObs <- as.numeric(is.na(y[,1]) | is.na(y[,2]))
  y[noObs==1,] = 0
  
  checkFit <- 1
  thetaAttempt <- theta
  while(attempts > 0 & checkFit) {
    if (!is.null(initialSANN) & method!='SANN') {
      #browser()
      message("Beginning SANN initialization ...")
      init <- optim(thetaAttempt, crwN2ll, method='SANN', control=initialSANN,
                    fixPar=fixPar, y=y, noObs=noObs,
                    delta=c(diff(data$TimeNum), 1), #a=initial.state$a, P=initial.state$P,
                    mov.mf=mov.mf, err.mfX=err.mfX, err.mfY=err.mfY, rho=rho, activity=activity,
                    n.mov=n.mov, n.errX=n.errX, n.errY=n.errY,
                    driftMod=driftMod, prior=prior, need.hess=FALSE, constr=constr)
      #thetaAttempt <- init$par
    } else init <- list(par=thetaAttempt)
    #if(any(init$par<lower)) init$par[init$par<lower] <- lower[init$par<lower] + 0.000001
    #if(any(init$par>upper)) init$par[init$par>upper] <- upper[init$par>upper] - 0.000001
    message("Beginning likelihood optimization ...")
    mle <- try(optim(init$par, crwN2ll, method=method, hessian=need.hess,
                     lower=constr$lower, upper=constr$upper, control=control,					  
                     fixPar=fixPar, y=y, noObs=noObs,
                     delta=c(diff(data$TimeNum), 1), #a=initial.state$a, P=initial.state$P,
                     mov.mf=mov.mf, err.mfX=err.mfX, err.mfY=err.mfY, activity=activity,
                     n.mov=n.mov, n.errX=n.errX, n.errY=n.errY, rho=rho,
                     driftMod=driftMod, prior=prior, need.hess=need.hess, constr=constr), silent=TRUE)
    attempts <- attempts - 1
    checkFit = inherits(mle, 'try-error')
  }
  if(inherits(mle, 'try-error')) return(mle)
  else {
    par <- fixPar
    par[is.na(fixPar)] <- mle$par
    Cmat <- matrix(NA, n.par, n.par)
    C.tmp <- try(2 * solve(mle$hessian), silent=TRUE)
    if (inherits(C.tmp, "try-error")) {
      cat("\nCannot calculate covariance matrix\n\n")
    } else Cmat[is.na(fixPar), is.na(fixPar)] <- C.tmp
    se <- sqrt(diag(Cmat))
    ci.l <- par - 1.96 * se
    ci.u <- par + 1.96 * se
    
    out <- list(par=par, estPar=mle$par, se=se, ci=cbind(L=ci.l, U=ci.u), Cmat=Cmat,
                loglik=-mle$value / 2, aic=mle$value + 2 * sum(is.na(fixPar)),
                coord=coord,fixPar=fixPar,
                convergence=mle$convergence, message=mle$message,
                activity=activity, random.drift=drift,
                mov.model=mov.model, err.model=err.model, n.par=n.par, nms=nms,
                n.mov=n.mov, n.errX=n.errX, n.errY=n.errY,
                mov.mf=mov.mf, err.mfX=err.mfX, err.mfY=err.mfY, rho=rho,
                Time.name=Time.name, init=init, data=data,
                lower=constr$lower, upper=constr$upper, prior=prior, need.hess=need.hess,
                runTime=difftime(Sys.time(), st))
    attr(out,"epsg") <- p4$epsg
    attr(out,"proj4") <- p4$proj4string
    attr(out, "time.scale") = ts
    class(out) <- c("crwFit")
    return(out)
  }
}
