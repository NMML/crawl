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
