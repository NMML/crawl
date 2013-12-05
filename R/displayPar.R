#'Display the order of parameters along with fixed values and starting values
#'
#'This function takes the model spesification arguments to the \code{\link{crwMLE}} function and displays a table
#'with the parameter names in the order that \code{crwMLE} will use during model fitting. This is useful for specifying 
#'values for the \code{fixPar} or \code{theta} (starting values for free parameters) arguments. 
#'
#'@param mov.model formula object specifying the time indexed covariates for
#'movement parameters.
#'@param err.model A 2-element list of formula objects specifying the time
#'indexed covariates for location error parameters.
#'@param stop.model formula object giving the covariate for the stopping
#'portion of the model.
#'@param drift.model logical indicating whether or not to include a random
#'drift component.
#'@param data data.frame object containg telemetry and covariate data. A
#'\code{\link{SpatialPointsDataFrame}} object from the package 'sp' will also be accepted.
#'@param theta starting values for parameter optimization.
#'@param fixPar Values of parameters which are held fixed to the given value.
#' 
#'@return A data frame with the following columns
#' 
#'\item{ParNames}{The names of the parameters specified by the arguments.}
#'
#'\item{fixPar}{The values specified by the \code{fixPar} argument for fixed values of the parameters. In model fitting, 
#'these values will remain fixed and will not be estimated.}
#'
#'\item{thetaIndex}{This column provides the index of each element of the theta argument and to which parameter it corresponds.}
#'
#'\item{thetaStart}{If a value is given for the \code{theta} argument it will be placed in this column and its elements will 
#'correspond to the \code{thetaIdx} column.}
#' 
#'@author Devin S. Johnson
#'@seealso \code{demo(northernFurSealDemo)} for example.
#'  
#'@export 

displayPar <- function(mov.model=~1, err.model=NULL, stop.model=NULL, drift.model=FALSE, data, theta, fixPar){
  ## SET UP MODEL MATRICES AND PARAMETERS ##
  if(inherits(data, "SpatialPoints")) {  
    polar.coord <- "+proj=longlat" %in% strsplit(proj4string(data), " ")[[1]]	
    coordVals <- as.data.frame(coordinates(data))	
    coord <- names(coordVals)	
    data <- cbind(slot(data,"data"), coordVals)    
  }
  errMod <- !is.null(err.model)
  stopMod <- !is.null(stop.model)
  driftMod <- drift.model
  mov.mf <- model.matrix(mov.model, model.frame(mov.model, data, na.action=na.pass))
  if (any(is.na(mov.mf))) stop("\nMissing values are not allowed in movement covariates!\n")
  n.mov <- ncol(mov.mf)
  if (errMod) {
    if (length(err.model) > 1) {
      err.mfY <- model.matrix(err.model[[2]],
                              model.frame(err.model[[2]], data, na.action=na.pass))
      err.mfY <- ifelse(is.na(err.mfY), 0, err.mfY)
      n.errY <- ncol(err.mfY)
    } else {
      err.mfY <- NULL
      n.errY <- 0
    }
    err.mfX <- model.matrix(err.model[[1]],
                            model.frame(err.model[[1]], data, na.action=na.pass))
    err.mfX <- ifelse(is.na(err.mfX), 0, err.mfX)
    n.errX <- ncol(err.mfX)
  } else {
    n.errY <- n.errX <- 0
    err.mfX <- err.mfY <- NULL
  }
  if (stopMod) {
    #stop.model
    stop.mf <- model.matrix(stop.model,
                            model.frame(stop.model, data, na.action=na.pass))
    if (ncol(stop.mf) > 2) stop("\nThere can only be one stopping variable >0 and <1\n")
    stop.mf <- as.double(stop.mf[, 2])
    if (any(stop.mf < 0) | any(stop.mf > 1)) stop("\nStop variable must be >0 and <1\n")
    if (any(is.na(stop.mf))) stop("\nMissing values are not allowed in the stopping variable!\n")
    n.stop <- 1
  } else stop.mf <- NULL
  n.drift <- as.integer(driftMod)
  n.stop <- as.integer(stopMod)
  b.nms <- paste("ln beta ", colnames(mov.mf), sep="")
  sig.nms <- paste("ln sigma ", colnames(mov.mf), sep="")
  if (errMod) {
    if (length(err.model) > 1) {
      tau.nms <- c(paste("ln tau.x ", colnames(err.mfX), sep=""),
                   paste("ln tau.y ", colnames(err.mfY), sep=""))
    } else tau.nms <- paste("ln tau ", colnames(err.mfX), sep="")
  } else tau.nms <- NULL
  if (stopMod) {stop.nms <- "ln phi"} else stop.nms <- NULL
  if (driftMod) {
    drift.nms <- c("ln sigma.drift/sigma", "ln psi-1")
  } else drift.nms <- NULL
  nms <- c(tau.nms, sig.nms, b.nms, stop.nms, drift.nms)
  n.par <- length(nms)
  if (missing(fixPar)) fixPar <- rep(NA, n.par)
  if (length(fixPar)!=n.par) stop("'fixPar' argument is not the right length! The number of parameters in the model is ", n.par, "\n")
  if (!missing(theta)) if(length(theta) != sum(is.na(fixPar))) stop("\nWrong number of parameters specified in start value.\n")
  thetaIdx <- fixPar
  thetaIdx[is.na(fixPar)] <- 1:sum(is.na(fixPar))
  thetaIdx[!is.na(fixPar)] <- NA
  out <- data.frame(ParNames=nms, fixPar=fixPar, thetaIdx=thetaIdx)
  if(!missing(theta)){
    thetaStart <- thetaIdx
    thetaStart[!is.na(thetaStart)] <- theta
    out <- cbind(out, thetaStart=thetaStart)
  }
  return(out)
}