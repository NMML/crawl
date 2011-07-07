"crwPredict" <- function(object.crwFit, predTime=NULL, speedEst=FALSE, flat=TRUE)
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
	if(inherits(predTime, "POSIXct")) predTime <- as.numeric(predTime)/3600

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
        if (stopMod) stop.mf <- as.matrix(expandPred(x=stop.mf, Time=origTime, predTime=predTime))
        if (!is.null(err.mfX)) err.mfX <- as.matrix(expandPred(x=err.mfX, Time=origTime, predTime=predTime))
        if (!is.null(err.mfY)) err.mfY <- as.matrix(expandPred(x=err.mfY, Time=origTime, predTime=predTime))
    }
    if (object.crwFit$polar.coord) {
        lonAdjVals <- cos(round(approx(data[, tn], data[, object.crwFit$coord[2]],
                                       data[, tn])$y, 0) * pi / 180)
        if(is.na(lonAdjVals[1])) stop("Error in Longitude correction: Check to see that prediction times do not predate first observation")
        if(is.na(lonAdjVals[nrow(data)])) {
         warning('Forcasting locations past last observation: Longitude correction will be based on the last observation')
          lonAdjVals <- lonAdjVals[!is.na(lonAdjVals)][cumsum(!is.na(lonAdjVals))]
        }
    } else lonAdjVals <- rep(1, nrow(data))
    data$locType[data[,tn]%in%predTime] <- 'p'
    delta <- c(diff(data[, tn]), 1)
    a1.x <- object.crwFit$initial.state$a1.x
    P1.x <- object.crwFit$initial.state$P1.x
    a1.y <- object.crwFit$initial.state$a1.y
    P1.y <- object.crwFit$initial.state$P1.y
    y <- data[, object.crwFit$coord[2]]
    x <- data[, object.crwFit$coord[1]]
    loctype <- ifelse(is.na(x) | is.na(y), 1, 0)
    y <- ifelse(loctype == 1, 9999, y)
    x <- ifelse(loctype == 1, 9999, x)
    N <- length(loctype)
   ###
   ### Process parameters for Fortran
   ###
   if (!is.null(err.mfX)) {
      theta.errX <- par[1:n.errX]
      tau2x <- exp(2 * err.mfX %*% theta.errX)
   } else tau2x <- rep(0.0, N)
   if (!is.null(err.mfY)) {
      theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
      tau2y <- exp(2 * err.mfY %*% theta.errY)
   } else tau2y <- tau2x
   theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
   sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
   b <- exp(mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)])
   stay <- rep(0, N)
   if (stopMod) {
      stop.mf <- stop.mf
      theta.stop <- par[(n.errX + n.errY + 2 * n.mov + 1)]
      b <- b / ((1 - stop.mf) ^ exp(theta.stop))
      stay <- ifelse(b==Inf, 1, 0)
      b <- ifelse(b==Inf, 9999, b) 
   }
   if (driftMod) {
      theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                                    (n.errX + n.errY + 2 * n.mov + 2)]
      b.drift <- exp(log(b) - log(1+exp(theta.drift[2])))
      sig2.drift <- exp(log(sig2) + 2 * theta.drift[1])
      call.predict <- "crwdrift_predict"
   } else {
      b.drift <- sig2.drift <- 0.0
      call.predict <- "crw_predict"
   }
                    
   movMats <- getQT(sig2, b, sig2.drift, b.drift, delta, driftMod)
   
    out <- .Fortran(call.predict,
                    tau2y=as.double(tau2y),
                    tau2x=as.double(tau2x),
                    Qmat=as.double(movMats$Qmat),
                    Tmat=as.double(movMats$Tmat),
                    x=as.double(x),
                    y=as.double(y),
                    loctype=as.integer(loctype),
                    stay=as.integer(stay),
                    ay=as.double(a1.y),
                    ax=as.double(a1.x),
                    Py=as.double(P1.y),
                    Px=as.double(P1.x),
                    lonadj=as.double(lonAdjVals),
                    N=as.integer(N),
                    lly=as.double(0),
                    llx=as.double(0),
                    My=as.double(rep(0,N)),
                    uy=as.double(rep(0,N)),
                    jky=as.double(rep(0,N)),
                    Mx=as.double(rep(0,N)),
                    ux=as.double(rep(0,N)),
                    jkx=as.double(rep(0,N)),
                    predy=as.double(matrix(0, N, 2 + driftMod)),
                    predx=as.double(matrix(0, N, 2 + driftMod)),
                    vary=as.double(array(0, c(2 + driftMod, 2 + driftMod, N))),
                    varx=as.double(array(0, c(2 + driftMod, 2 + driftMod, N))),
                    package="crawl")
                    
    predy <- data.frame(matrix(out$predy, N, 2 + driftMod))
    if (driftMod) {
        names(predy) <- c("mu.y", "theta.y", "gamma.y")
    } else names(predy) <- c("mu.y", "nu.y")
    predx <- data.frame(matrix(out$predx, N, 2 + driftMod))
    if (driftMod) {
        names(predx) <- c("mu.x", "theta.x", "gamma.x")
    } else names(predx) <- c("mu.x", "nu.x")
    vary <- zapsmall(array(out$vary, c(2 + driftMod, 2 + driftMod, N)))
    varx <- zapsmall(array(out$varx, c(2 + driftMod, 2 + driftMod, N)))
    obsFit <- data.frame(shock.x=ifelse(out$Mx!=0,out$ux/out$Mx,NA), 
                         predObs.x=ifelse(out$Mx!=0,out$jkx,NA), 
                         Vshock.x=ifelse(out$Mx!=0,1/out$Mx,NA), 
                         shock.y=ifelse(out$My!=0,out$uy/out$My,NA), 
                         predObs.y=ifelse(out$My!=0,out$jky,NA), 
                         Vshock.y=ifelse(out$My!=0,1/out$My,NA))
    obsFit$outlier.chisq <- out$ux^2/out$Mx + out$uy^2/out$My
    obsFit$naive.p.val <- 1 - pchisq(obsFit$outlier.chisq, 2)
    if (speedEst) {
        log.speed <- logSpeed(predx, predy, varx, vary, object.crwFit$polar.coord)
    } else log.speed <- NULL
    out <- list(originalData=fillCols(data), alpha.hat.y=predy, alpha.hat.x=predx,
                V.hat.y=vary, V.hat.x=varx, speed=log.speed, loglik=out$lly + out$llx)
    if (flat) {
        out <- cbind(fillCols(as.flat(out)), obsFit)
        attr(out, "flat") <- TRUE
        attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    	attr(out, "random.drift") <- driftMod
    	attr(out, "stop.model") <- object.crwFit$stopMod
    	attr(out, "polar.coord") <- object.crwFit$polar.coord
    	attr(out, "Time.name") <- tn
    } else {
       out <- append(out, list(fit.test=obsFit))
        attr(out, "flat") <- FALSE
        attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    	attr(out, "random.drift") <- driftMod
    	attr(out, "stop.model") <- object.crwFit$stopMod
    	attr(out, "polar.coord") <- object.crwFit$polar.coord
    	attr(out, "Time.name") <- tn
    }
    class(out) <- c(class(out),"crwPredict")
    return(out)
}
