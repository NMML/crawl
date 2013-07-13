#' -2 * log-likelihood for CTCRW models
#' 
#' This function is designed for primary use within the \code{\link{crwMLE}}
#' model fitting function. But, it can be accessed for advanced \code{R} and
#' \code{crawl} users. Uses the state-space parameterization and Kalman filter
#' method presented in Johnson et al. (2008).
#' 

#' 
#' This function calls compiled Fortran code which can be viewed in the
#' \code{src} directory of the crawl library.
#' 
#' @param theta parameter values.
#' @param fixPar values of parameters held fixed (contains \code{NA} for
#' \code{theta} values).
#' @param y latitude locations.
#' @param x longitude loations.
#' @param loctype vector with 1 for unobserved locations, and 0 for observed locations.
#' @param delta time difference to next location.
#' @param a1.y initial state value for latitude.
#' @param a1.x initial state value for longitude.
#' @param P1.x intial state covariance matrix for latitude.
#' @param P1.y inital state covariance matrix for longitude.
#' @param lonAdj = 1/cos(y*pi/180) as described in Johnson et al. (2008) for
#' polar coords. = 1 for non-polar coords.
#' @param mov.mf Movement covariate data.
#' @param err.mfX longitude error covariate data.
#' @param err.mfY latitude error covariate data.
#' @param stop.mf stopping covariate.
#' @param n.errX number or longitude error parameters.
#' @param n.errY number of latitude error parameters.
#' @param n.mov number or movement parameters.
#' @param stopMod Logical. indicates wheteher a stop model is specified.
#' @param driftMod Logical. inicates whether a drift model is specified.
#' @param prior Function of theta that returns the log-density of the prior
#' @param need.hess Whether or not the Hessian will need to be calculated from
#' this call
#' @param constr Named list giving the parameter constraints
#' @return -2 * log-likelihood value for specified CTCRW model.
#' @author Devin S. Johnson
#' @seealso \code{\link{crwMLE}}
#' @references Johnson, D., J. London, M. -A. Lea, and J. Durban. 2008.
#' Continuous-time model for animal telemetry data. Ecology 89:1208-1215.
#' @export

"crwN2ll" <- function(theta, fixPar, y, x, loctype, delta, a1.y, a1.x,
                      P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
                      n.errX, n.errY, n.mov, stopMod, driftMod, prior, need.hess, constr=list(lower=-Inf, upper=Inf))
{
	if(!need.hess & any(theta < constr$lower | theta > constr$upper)) return(Inf)
    N <- length(y)
    par <- fixPar
    par[is.na(fixPar)] <- theta
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
      sig2.drift <- exp(log(sig2) + 2 * theta.drift[1]) #rep(exp(2*theta.drift[1]), length(sig2)) #
      call.lik <- "crwdriftn2ll"
   } else {
      b.drift <- sig2.drift <- 0.0
      call.lik <- "crwn2ll"
   }
   movMats <- getQT(sig2, b, sig2.drift, b.drift, delta, driftMod)
    out <- .Fortran(call.lik,
                    tau2y=as.double(tau2y),
                    tau2x=as.double(tau2x),
                    Qmat=as.double(as.matrix(movMats$Qmat)),
                    Tmat=as.double(as.matrix(movMats$Tmat)),
                    x=as.double(x),
                    y=as.double(y),
                    loctype=as.integer(loctype),
                    stay=as.integer(stay),
                    ay=as.double(a1.y),
                    ax=as.double(a1.x),
                    Py=as.double(P1.y),
                    Px=as.double(P1.x),
                    lonadj=as.double(lonAdj),
                    N=as.integer(N),
                    lly=as.double(0),
                    llx=as.double(0),
                    package="crawl")
    if(is.null(prior)) return(-2 * (out$lly + out$llx))
	else return(-2 * (out$lly + out$llx + prior(theta)))
}
