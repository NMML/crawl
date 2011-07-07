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