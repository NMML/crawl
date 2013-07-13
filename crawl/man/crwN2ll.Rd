\name{crwN2ll}
\alias{crwN2ll}
\title{-2 * log-likelihood for CTCRW models}
\usage{
  crwN2ll(theta, fixPar, y, x, loctype, delta, a1.y, a1.x,
    P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
    n.errX, n.errY, n.mov, stopMod, driftMod, prior,
    need.hess, constr = list(lower = -Inf, upper = Inf))
}
\arguments{
  \item{theta}{parameter values.}

  \item{fixPar}{values of parameters held fixed (contains
  \code{NA} for \code{theta} values).}

  \item{y}{latitude locations.}

  \item{x}{longitude loations.}

  \item{loctype}{vector with 1 for unobserved locations,
  and 0 for observed locations.}

  \item{delta}{time difference to next location.}

  \item{a1.y}{initial state value for latitude.}

  \item{a1.x}{initial state value for longitude.}

  \item{P1.x}{intial state covariance matrix for latitude.}

  \item{P1.y}{inital state covariance matrix for
  longitude.}

  \item{lonAdj}{= 1/cos(y*pi/180) as described in Johnson
  et al. (2008) for polar coords. = 1 for non-polar
  coords.}

  \item{mov.mf}{Movement covariate data.}

  \item{err.mfX}{longitude error covariate data.}

  \item{err.mfY}{latitude error covariate data.}

  \item{stop.mf}{stopping covariate.}

  \item{n.errX}{number or longitude error parameters.}

  \item{n.errY}{number of latitude error parameters.}

  \item{n.mov}{number or movement parameters.}

  \item{stopMod}{Logical. indicates wheteher a stop model
  is specified.}

  \item{driftMod}{Logical. inicates whether a drift model
  is specified.}

  \item{prior}{Function of theta that returns the
  log-density of the prior}

  \item{need.hess}{Whether or not the Hessian will need to
  be calculated from this call}

  \item{constr}{Named list giving the parameter
  constraints}
}
\value{
  -2 * log-likelihood value for specified CTCRW model.
}
\description{
  This function is designed for primary use within the
  \code{\link{crwMLE}} model fitting function. But, it can
  be accessed for advanced \code{R} and \code{crawl} users.
  Uses the state-space parameterization and Kalman filter
  method presented in Johnson et al. (2008).
}
\details{
  This function calls compiled Fortran code which can be
  viewed in the \code{src} directory of the crawl library.
}
\author{
  Devin S. Johnson
}
\references{
  Johnson, D., J. London, M. -A. Lea, and J. Durban. 2008.
  Continuous-time model for animal telemetry data. Ecology
  89:1208-1215.
}
\seealso{
  \code{\link{crwMLE}}
}

