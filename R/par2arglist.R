
#' @author Devin Johnson
#' @importFrom stats var
#' @importFrom utils tail

par2arglist = function(theta, fixPar, y, noObs, delta, 
                       mov.mf, err.mfX, err.mfY, rho, activity,
                       n.errX, n.errY, n.mov, driftMod){
  N <- nrow(y)
  par <- fixPar
  par[is.na(fixPar)] <- theta
  
  out = vector("list",8)
  names(out) = c("Hmat", "b", "sig2", "active", "b.drift", "sig2.drift", "a", "P")

  # Hmat
  if (!is.null(err.mfX)) {
    theta.errX <- par[1:n.errX]
    Hmat <- exp(2 * err.mfX %*% theta.errX)
  } else Hmat <- rep(0.0, N)
  if (!is.null(err.mfY)) {
    theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
    Hmat <- cbind(Hmat,exp(2 * err.mfY %*% theta.errY))
  } else Hmat <- cbind(Hmat, Hmat)
  if(!is.null(rho)){
    Hmat = cbind(Hmat, exp(log(Hmat[,1])/2 + log(Hmat[,2])/2)*rho)
  } else {Hmat = cbind(Hmat, rep(0,N))}
  Hmat[noObs==1,] = 0
  out$Hmat = Hmat
  
  # b, sig2, active
  theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
  out$sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
  out$b <- exp(mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)])
  if (!is.null(activity)) {
    theta.stop <- par[(n.errX + n.errY + 2 * n.mov + 1)]
    out$b <- out$b / ((activity) ^ exp(theta.stop))
    out$sig2 = out$sig2 * ((activity) ^ exp(theta.stop))
    out$active <- ifelse(out$b==Inf, 0, 1)
    out$b <- ifelse(out$b==Inf, 0, out$b) 
  } else {
    out$active=rep(1,N)
  }
  
  # b.drift, sig2.drift, a, and P
  if (driftMod) {
    theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                         (n.errX + n.errY + 2 * n.mov + 2)]
    out$b.drift <- exp(log(out$b) - log(1+exp(theta.drift[2])))
    out$sig2.drift <- exp(log(out$sig2) + 2 * theta.drift[1]) 
    out$a = c(y[1,1], 0, 0, y[1,2],0, 0)
    out$P = diag(c(var(y[noObs==0,1], na.rm=T), out$sig2[1]*out$b[1], out$sig2.drift[1]*out$b.drift[1], 
               var(y[noObs==0,2],na.rm=T), out$sig2[1]*out$b[1], out$sig2.drift[1]*out$b.drift[1]))
  } else {
    out$b.drift = NULL
    out$sig2.drift = NULL
    out$a = c(y[1,1], 0, y[1,2],0)
    out$P = diag(c(var(y[noObs==0,1], na.rm=T), out$sig2[1]*out$b[1], var(y[noObs==0,2],na.rm=T), out$sig2[1]*out$b[1]))
  }
  return(out)
}