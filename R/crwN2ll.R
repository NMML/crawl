#' -2 * log-likelihood for CTCRW models
#' 
#' This function is designed for primary use within the \code{\link{crwMLE}}
#' model fitting function. But, it can be accessed for advanced \code{R} and
#' \code{crawl} users. Uses the state-space parameterization and Kalman filter
#' method presented in Johnson et al. (2008).
#' 

#' 
#' This function calls compiled C++ code which can be viewed in the
#' \code{src} directory of the crawl source package.
#' 
#' @param theta parameter values.
#' @param fixPar values of parameters held fixed (contains \code{NA} for
#' \code{theta} values).
#' @param y N by 2 matrix of coordinates with the longitude coordinate in the first column.
#' @param noObs vector with 1 for unobserved locations, and 0 for observed locations.
#' @param delta time difference to next location.
#' @param mov.mf Movement covariate data.
#' @param err.mfX longitude error covariate data.
#' @param err.mfY latitude error covariate data.
#' @param rho A vector of known correlation coefficients for the error model, typically used for modern ARGOS data.
#' @param activity Stopping covariate (= 0 if animal is not moving).
#' @param n.errX number or longitude error parameters.
#' @param n.errY number of latitude error parameters.
#' @param n.mov number or movement parameters.
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

crwN2ll = function(theta, fixPar, y, noObs, delta, #a, P, 
                   mov.mf, err.mfX, err.mfY, rho=NULL, activity=NULL,
                   n.errX, n.errY, n.mov, driftMod, prior, need.hess, 
                   constr=list(lower=-Inf, upper=Inf))
{
  if(!need.hess & any(theta < constr$lower | theta > constr$upper)) return(Inf)
  N <- nrow(y)
  par <- fixPar
  par[is.na(fixPar)] <- theta
  
  argslist = par2arglist(theta, fixPar, y, noObs, delta,
                          mov.mf, err.mfX, err.mfY, rho=NULL, activity=NULL,
                          n.errX, n.errY, n.mov, driftMod)
  
  if (driftMod) {
    ll <- CTCRWNLL_DRIFT(as.matrix(y), argslist$Hmat, argslist$b, argslist$b.drift, 
                          argslist$sig2, argslist$sig2.drift, delta, noObs, argslist$active, argslist$a,  argslist$P)$ll
  } else {
    ll <- CTCRWNLL(as.matrix(y), argslist$Hmat, argslist$b, argslist$sig2, delta, noObs, argslist$active, argslist$a,  argslist$P)$ll
  }

  if(is.null(prior)){
    ll = -2 * ll
  } else {
    ll = -2 * (ll + prior(theta))
  }
  return(ll)
}
