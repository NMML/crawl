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
#' A full model specification involves 4 components: a movement model, a
#' stopping model, 2 location error models, and a drift indication. The
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
#' process, a drift variance and a beta multiplier. If \code{polar.coord=TRUE}
#' then the ad-hoc logitude correction factor described by Johnson et al.
#' (2008) (Ecology 89:1208-1215) is used to adjust the variance scale for the
#' longitude mdoel.
#' 
#' The \code{inital.state} is a list with the following elemets (with the exact
#' names):
#' 
#' \code{a1.y} A vector with initial state values for the \dQuote{latitude}
#' coordinate. It have 2 elemets (location at time 1, velocity at time 1) for
#' non-drift models and 3 elemets for drift models (location at time 1,
#' velocity at time 1, drift velocity at time 1) for driftmodels,
#' 
#' \code{P1.y} Covarince matrix for the state at time 1 (measure of uncertainty
#' for your inital state) \code{a1.y},
#' 
#' \code{a1.x} Same as \code{a1.y}, but in the \dQuote{longitude} coordinate,
#' 
#' \code{P1.x} Same as \code{P1.y}, but in the \dQuote{longitude} coordinate.
#' 
#' \code{theta} and \code{fixPar} are vectors with the appropriate number or
#' parameters. \code{theta} contains only those paraemters which are to be
#' estimated, while \code{fixPar} contains all parameter values with \code{NA}
#' for parameters which are to be estimated.
#' 
#' The data set specified by \code{data} must contain a numeric or POSIXct
#' column which is used as the time index for analysis. The column name is
#' specified by the \code{Time.name} argument. If a POSIXct column is used it
#' is internally converted to a numeric vector with units of hours. If this is
#' not appropriate for your data, it is better to convert it yourself prior to
#' analysis with crawl. Also, for stopping models, the stopping covariate must
#' be between 0 and 1 inclusive, with 1 representing complete stop of the
#' animal (no true movement, however, location error can still occur) and 0
#' represent unhindered movement. The coordinate location should have \code{NA}
#' where no location is recorded, but there is a change in the movment
#' covariates.
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
#' @param stop.model formula object giving the covariate for the stopping
#' portion of the model.
#' @param drift.model logical indicating whether or not to include a random
#' drift component.
#' @param data data.frame object containg telemetry and covariate data. A
#' 'SpatialPointsDataFrame' object from the package 'sp' will also be accepted.
#' In which case the \code{polar.coord} and \code{coord} values will be taken
#' from the spatial data set and ignored in the arguments.
#' @param coord A 2-vector of character values giving the names of the "X" and
#' "Y" coordinates in \code{data}.
#' @param polar.coord logical indicating location are in degrees latitude and
#' longitude.
#' @param Time.name character indicating name of the location time column
#' @param initial.state list object containg the inital state of the Kalman
#' filter.
#' @param theta starting values for parameter optimization.
#' @param fixPar Values of parameters which are held fixed to the given value.
#' @param method Optimization method that is passed to \code{\link{optim}}.
#' @param control Control list which is passed to \code{\link{optim}}.
#' @param constr Named list with elements \code{lower} and \code{upper} that
#' are vectors the same length as theta giving the box constraints for the
#' parameters
#' @param prior A function returning the log-density function of the parameter
#' prior specification
#' @param need.hess A logical value which decides whether or not to evaluate
#' the Hessian for parameter standard errors
#' @param initialSANN Control list for \code{\link{optim}} when simulated
#' annealing is used for obtaining start values. See details
#' @param attempts The number of times likelihood optimization will be
#' attempted
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
#' \item{initial.state}{Intial state provided to \code{crwMLE} for model
#' fitting}
#' 
#' \item{coord}{Coordinate names provided for fitting}
#' 
#' \item{fixPar}{Fixed parameter values provided}
#' 
#' \item{convergence}{Indicator of convergence (0 = converged)}
#' 
#' \item{message}{Meesages given by \code{optim} during parameter optimization}
#' 
#' \item{stop.model}{Model provided for stopping variable}
#' 
#' \item{random.drift}{Logical value indicating random drift model}
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
#' @author Devin S. Johnson
#' @seealso \code{\link{northernFurSeal}} for additional examples.
#' @examples
#' 
#' 
#' data(harborSeal)
#' head(harborSeal)
#' ## Calculate Log multipliers for Argos error
#' argosClasses <- c("3", "2", "1", "0", "A", "B")
#' ArgosMultFactors <- data.frame(Argos_loc_class=argosClasses,
#'                                errX=log(c(1, 1.5, 4, 14, 5.21, 20.78)),
#'                                errY=log(c(1, 1.5, 4, 14, 11.08, 31.03)))
#' hsNew <- merge(harborSeal, ArgosMultFactors, by=c("Argos_loc_class"), all=TRUE)
#' hsNew <- hsNew[order(hsNew$Time), ]
#' head(hsNew)
#' 
#' ## Initial state values
#' initial.dry <- list(
#'   a1.x=c(harborSeal$longitude[1],0),
#'   a1.y=c(harborSeal$latitude[1],0),
#'   P1.x=diag(c(1,1)),
#'   P1.y=diag(c(1,1))
#' )
#' 
#' ##Fit model as given in Johnson et al. (2008) Ecology 89:1208-1215
#' ## Start values for theta come from the estimates in Johnson et al. (2008)
#' 
#' fit1 <- crwMLE(
#'   mov.model=~1, err.model=list(x=~errX, y=~errY), stop.model=~DryTime,
#'   data=hsNew, coord=c("longitude","latitude"), polar.coord=TRUE, Time.name="Time", 
#'   initial.state=initial.dry, fixPar=c(NA, 1, NA, 1, NA, NA, NA), theta=c(-6,-7,-4,-0.5,-1),
#'   control=list(maxit=2000, trace=1, REPORT=1),
#' )
#' 
#' fit1
#' str(fit1)
#' 
#' ##Use simulated annealing to obtain start values and place constraints on the parameters
#' 
#' set.seed(123)
#' fit2 <- crwMLE(
#'   mov.model=~1, err.model=list(x=~errX, y=~errY), stop.model=~DryTime,
#'   data=hsNew, coord=c("longitude","latitude"), Time.name="Time", polar.coord=TRUE,
#'   initial.state=initial.dry, fixPar=c(NA, 1, NA, 1, NA, NA, NA),
#'   control=list(maxit=2000, trace=1, REPORT=1),
#'   constr=list(lower=c(-6, -6, -Inf, -Inf, -Inf), upper=Inf),
#'   initialSANN=list(maxit=100, temp=5, tmax=5, trace=1, REPORT=2)
#' )
#' 
#' fit2
#' 
#' ##See simulated annealing start values
#' fit2$init$par
#' @export
"crwMLE" <- function(mov.model=~1, err.model=NULL, stop.model=NULL, drift.model=FALSE,
                     data, coord=c("x", "y"), polar.coord, Time.name,
                     initial.state, theta, fixPar, method="L-BFGS-B", control=NULL, constr=list(lower=-Inf, upper=Inf), 
					 prior=NULL, need.hess=TRUE, initialSANN=NULL, attempts=1)
{
    st <- Sys.time()
    if (missing(Time.name)) stop("Argument 'Time.name' missing. Please specify")

    
    ### Transform 'sp' package SpatialPointsDataFrame
    if(inherits(data, "trip")){
		Time.name <- data@TOR.columns[1]
	}
    if(inherits(data, "SpatialPoints")) {	
    	polar.coord <- "+proj=longlat" %in% strsplit(proj4string(data), " ")[[1]]	
    	coordVals <- as.data.frame(coordinates(data))	
    	coord <- names(coordVals)	
    	data <- cbind(slot(data,"data"), coordVals)    
    }
	if(inherits(data[,Time.name],"POSIXct")){
		data$TimeNum <- as.numeric(data[,Time.name])/3600
		Time.name <- "TimeNum"
	}
	
	
	### Check for duplicate time records ###
	if(any(diff(data[,Time.name])==0)) stop("\nERROR: There are duplicate time records for some data entries! Please remove before proceeding.\n")
	
	
    ## SET UP MODEL MATRICES AND PARAMETERS ##
    errMod <- !is.null(err.model)
    stopMod <- !is.null(stop.model)
    driftMod <- drift.model
    if (
        length(initial.state$a1.y) != driftMod+2 |
        length(initial.state$a1.x) != driftMod+2 |
        all(dim(initial.state$P1.y) != c(driftMod+2, driftMod+2)) |
        all(dim(initial.state$P1.x) != c(driftMod+2, driftMod+2))
        ) stop("Dimentions of 'initial.state' argument are not correct for the specified model")
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
	if(!(length(constr$lower)==1 | length(constr$lower)==sum(is.na(fixPar)))) stop("The number of lower contraints specified is not correct! The number of free parameters is ", sum(is.na(fixPar)),"\n")
	if(!(length(constr$upper)==1 | length(constr$upper)==sum(is.na(fixPar)))) stop("The number of upper contraints specified is not correct! The number of free parameters is ", sum(is.na(fixPar)),"\n")
    if(length(constr$upper)==1) constr$upper <- rep(constr$upper, sum(is.na(fixPar)))
	if(length(constr$lower)==1) constr$lower <- rep(constr$lower, sum(is.na(fixPar)))
	if (missing(theta)) theta <- ifelse(constr$lower > -Inf, constr$lower+0.001, 0.0)
    #theta <- ifelse(is.na(theta), 0.00001, theta)
    if(driftMod & is.na(fixPar[n.par])) theta[sum(is.na(fixPar))] <- log(diff(range(data[,Time.name]))/9)
    if (length(theta) != sum(is.na(fixPar))) {
        stop("\nWrong number of parameters specified in start value.\n")
    }

    ## PROCESS DATA AND LONGITUDE ADJUSTMENT FOR POLAR COORDS ##
    x <- as.vector(data[, coord[1]])
    y <- as.vector(data[, coord[2]])
    loctype <- is.na(x) | is.na(y)
    y.lik <- ifelse(is.na(y), 9999, y)
    x.lik <- ifelse(is.na(x), 9999, x)
    if (polar.coord) {
        lonAdjVals <- cos(round(approx(data[, Time.name], data[, coord[2]],
                                       data[, Time.name])$y, 0) * pi / 180)
    } else lonAdjVals <- rep(1, nrow(data))

    ## DEFINING OPTIMIZATION PROCEDURE ##
#     if(missing(lower)) {
#     	if (method=='SANN' | !driftMod) lower <- -Inf
#     	else {
#       		if (is.na(fixPar[n.par])) {
#           		lower <- c(rep(-Inf, length(theta)-1), 0)
#       		} else lower <- -Inf 
#     	}
#     }

    checkFit <- 1
    thetaAttempt <- theta
    while(attempts > 0 & checkFit == 1) {
      if (!is.null(initialSANN) & method!='SANN') {
         init <- optim(thetaAttempt, crwN2ll, method='SANN', control=initialSANN,
                       fixPar=fixPar, y=y.lik, x=x.lik, loctype=loctype,
                       delta=c(diff(data[, Time.name]), 1), a1.y=initial.state$a1.y,
                       a1.x=initial.state$a1.x, P1.x=initial.state$P1.x,
                       P1.y=initial.state$P1.y, lonAdj=lonAdjVals, mov.mf=mov.mf,
                       err.mfX=err.mfX, err.mfY=err.mfY, stop.mf=stop.mf,
                       n.mov=n.mov, n.errX=n.errX, n.errY=n.errY, stopMod=stopMod,
                       driftMod=driftMod, prior=prior, need.hess=FALSE, constr=constr)
         #thetaAttempt <- init$par
      } else init <- list(par=thetaAttempt)
      #if(any(init$par<lower)) init$par[init$par<lower] <- lower[init$par<lower] + 0.000001
      #if(any(init$par>upper)) init$par[init$par>upper] <- upper[init$par>upper] - 0.000001
      mle <- try(optim(init$par, crwN2ll, method=method, hessian=need.hess,
				   lower=constr$lower, upper=constr$upper, control=control,					  
                   fixPar=fixPar, y=y.lik, x=x.lik, loctype=loctype,
                   delta=c(diff(data[, Time.name]), 1), a1.y=initial.state$a1.y,
                   a1.x=initial.state$a1.x, P1.x=initial.state$P1.x,
                   P1.y=initial.state$P1.y, lonAdj=lonAdjVals, mov.mf=mov.mf,
                   err.mfX=err.mfX, err.mfY=err.mfY, stop.mf=stop.mf,
                   n.mov=n.mov, n.errX=n.errX, n.errY=n.errY, stopMod=stopMod,
                   driftMod=driftMod, prior=prior, need.hess=need.hess), silent=TRUE)
      attempts <- attempts - 1
      checkFit <- 1.0*(inherits(mle, 'try-error'))
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
                   initial.state=initial.state, coord=coord, fixPar=fixPar,
                   convergence=mle$convergence, message=mle$message,
                   stop.model=stop.model, random.drift=drift.model,
                   mov.model=mov.model, err.model=err.model, n.par=n.par, nms=nms,
                   n.mov=n.mov, n.errX=n.errX, n.errY=n.errY,
                   mov.mf=mov.mf, err.mfX=err.mfX, err.mfY=err.mfY, stop.mf=stop.mf,
                   polar.coord=polar.coord, Time.name=Time.name, init=init, data=data,
                   lower=constr$lower, upper=constr$upper, prior=prior, need.hess=need.hess,
                   runTime=difftime(Sys.time(), st))
       class(out) <- c("crwFit")
       return(out)
    }
}
