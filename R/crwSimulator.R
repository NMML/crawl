#' Construct a posterior simulation object for the CTCRW state vectors
#' 

#' 
#' The \code{crwSimulator} function uses a fitted model object from
#' \code{crwMLE} and a set of prediction times to construct a list from which
#' \code{\link{crwPostIS}} will draw a sample from either the posterior
#' distribution of the state vectors conditional on fitted parameters or a full
#' posterior draw from an importance sample of the parameters.
#' 

#' 
#' The crwSimulator function produces a list and preprocesses the necessary
#' components for repeated track simulation from a fitted CTCRW model from
#' \code{\link{crwMLE}}. The \code{method} argument can be one of \code{"IS"}
#' or \code{"quadrature"}. If method="IS" is chosen standard importance
#' sampling will be used to calculate the appropriate weights via t proposal
#' with df degrees of freedom.  If df=Inf (default) then a multivariate normal
#' distribution is used to approximate the parameter posterior. If
#' \code{method="quadrature"}, then a regular grid over the posterior is used
#' to calculate the weights. The argument \code{grid.eps} controls the
#' quadrature grid. The arguments are approximately the upper and lower limit
#' in terms of standard deviations of the posterior.  The default is
#' \code{grid.eps}, in units of 1sd. If \code{object.crwFit} was fitted with
#' \code{crwArgoFilter}, then the returned list will also include \code{p.out},
#' which is the approximate probability that the observation is an outlier.
#' 
#' @param object.crwFit A model object from \code{\link{crwMLE}}.
#' @param predTime vector of additional prediction times.
#' @param method Method for obtaining weights for movement parameter samples
#' @param parIS Size of the parameter importance sample
#' @param df Degrees of freedom for the t approximation to the parameter
#' posterior
#' @param grid.eps Grid size for \code{method="quadrature"}
#' @param crit Criterion for deciding "significance" of quadrature points
#' (difference in log-likelihood)
#' @param scale Scale multiplier for the covariance matrix of the t
#' approximation
#' @return
#' 
#' List with the following elements:
#' 
#' \item{x}{Longitude coordinate with NA at prediction times}
#' 
#' \item{y}{Similar to above for latitude}
#' 
#' \item{locType}{Indicates prediction types with a "p" or observation times
#' with an "o"} \item{P1.y}{Initial state covariance for latitude}
#' 
#' \item{P1.x}{Initial state covariance for longitude}
#' 
#' \item{a1.y}{Initial latitude state}
#' 
#' \item{a1.x}{Initial longitude state}
#' 
#' \item{n.errX}{number of longitude error model parameters}
#' 
#' \item{n.errY}{number of latitude error model parameters}
#' 
#' \item{delta}{vector of time differences}
#' 
#' \item{driftMod}{Logical. indicates random drift model}
#' 
#' \item{stopMod}{Logical. Indicated stop model fitted}
#' 
#' \item{stop.mf}{stop model design matrix}
#' 
#' \item{err.mfX}{Longitude error model design matrix}
#' 
#' \item{err.mfY}{Latitude error model design matrix}
#' 
#' \item{mov.mf}{Movement model design matrix}
#' 
#' \item{fixPar}{Fixed values for parameters in model fitting}
#' 
#' \item{Cmat}{Covaraince matrix for parameter sampling distribution}
#' 
#' \item{Lmat}{Cholesky decomposition of Cmat}
#' 
#' \item{par}{fitted parameter values}
#' 
#' \item{N}{Total number of locations}
#' 
#' \item{loglik}{log likelihood of the fitted model}
#' 
#' \item{Time}{vector of observation times}
#' 
#' \item{coord}{names of coordinate vectors in original data}
#' 
#' \item{Time.name}{Name of the observation times vector in the original data}
#' 
#' \item{thetaSampList}{A list containing a data frame of parameter vectors and
#' their associated probabilities for a resample}
#' @author Devin S. Johnson
#' @seealso See \code{demo(northernFurSealDemo)} for example.
#' @export
`crwSimulator` <-
function(object.crwFit, predTime=NULL, method="IS", parIS=1000, df=Inf, grid.eps=1, crit=2.5, scale=1)
{
   ## Model definition/parameters ##
   data <- object.crwFit$data
   driftMod <- object.crwFit$random.drift
   stopMod <- !is.null(object.crwFit$stop.model)
   mov.mf <- object.crwFit$mov.mf
   stop.mf <- object.crwFit$stop.mf
   err.mfX <- object.crwFit$err.mfX
   err.mfY <- object.crwFit$err.mfY
   par <- object.crwFit$par
   n.errX <- object.crwFit$n.errX
   n.errY <- object.crwFit$n.errY
   n.mov <- object.crwFit$n.mov
   tn <- object.crwFit$Time.name
   if(any(diff(data[,tn])==0)) stop("\nERROR: There are duplicate time records for some data entries! Please remove before proceeding.\n")
   if(inherits(predTime, "POSIXct")){
	   predTime <- as.numeric(predTime)/3600
   }
   ## Data setup ##
   if (!is.null(predTime)) {
	  if(inherits(predTime, "POSIXct")) predTime <- as.numeric(predTime)/3600
      origTime <- data[, tn]
      data$predType <- "o"
      predData <- data.frame(predTime, "p")
      names(predData) <- c(tn, "predType")
      data <- merge(data, predData, by=c(tn, "predType"), all=TRUE)
      dups <- duplicated(data[, tn]) & data[,"predType"]=="p"
      data <- data[!dups, ]
      mov.mf <- as.matrix(expandPred(x=mov.mf, Time=origTime, predTime=predTime))
      if (stopMod) stop.mf <- as.matrix(expandPred(x=stop.mf, Time=origTime, predTime=predTime))
      if (!is.null(err.mfX)) err.mfX <- as.matrix(expandPred(x=err.mfX, Time=origTime, predTime=predTime))
      if (!is.null(err.mfY)) err.mfY <- as.matrix(expandPred(x=err.mfY, Time=origTime, predTime=predTime))
      data$predType[data[,tn]%in%predTime] <- 'p'
    }
    if (object.crwFit$polar.coord) {
        lonAdjVals <- cos(round(approx(data[, tn], data[, object.crwFit$coord[2]],
                                       data[, tn])$y, 0) * pi / 180)
    } else lonAdjVals <- rep(1, nrow(data))
    delta <- c(diff(data[, tn]), 1)
    a1.x <- object.crwFit$initial.state$a1.x
    P1.x <- object.crwFit$initial.state$P1.x
    a1.y <- object.crwFit$initial.state$a1.y
    P1.y <- object.crwFit$initial.state$P1.y
    y <- data[, object.crwFit$coord[2]]
    x <- data[, object.crwFit$coord[1]]
    if(is.null(data$predType)) data$predType <- ifelse(is.na(y)|is.na(x), "p", "o")
    loctype <- ifelse(is.na(x) | is.na(y), 1, 0)
    y <- ifelse(loctype == 1, 9999, y)
    x <- ifelse(loctype == 1, 9999, x)
    #Lmat <- t(chol(object.crwFit$Cmat[is.na(object.crwFit$fixPar),is.na(object.crwFit$fixPar)]))
    out <- list(x=x, y=y, loctype=loctype, P1.y=P1.y, P1.x=P1.x, a1.x=a1.x,
                a1.y=a1.y, n.errX=n.errX, n.errY=n.errY, n.mov=n.mov,
                delta=delta, lonAdj=lonAdjVals, driftMod=driftMod,
                stopMod=stopMod, stop.mf=stop.mf, err.mfX=err.mfX,
                err.mfY=err.mfY, mov.mf=mov.mf, fixPar=object.crwFit$fixPar,
                Cmat=object.crwFit$Cmat, predType=data$predType,
                par=object.crwFit$par, nms=object.crwFit$nms, N=nrow(data), lower=object.crwFit$lower, 
                upper=object.crwFit$upper, #lambda=lambdaOut,
                loglik=object.crwFit$loglik, data[,tn], 
                coord=object.crwFit$coord, Time.name=object.crwFit$Time.name, prior=object.crwFit$prior)
    names(out)[28] <- object.crwFit$Time.name
   class(out) <- 'crwSimulator'
   if(parIS>1 & object.crwFit$need.hess==TRUE) out <- crwSamplePar(out, method=method, size=parIS, df=df, grid.eps=grid.eps, crit=crit, scale=scale)
   return(out)
}
