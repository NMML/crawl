\name{logSpeed}
\alias{logSpeed}
\title{calculate estimate of log speed from a crwPredict object}
\arguments{
  \item{predx}{matrix of longitude state estimates.}

  \item{predy}{matrix of latitude state estimates.}

  \item{varx}{an array of covariance matrices for predx.}

  \item{vary}{ana array of covariance matrices for predy.}

  \item{polar.coord}{Logical. TRUE if coordinates are
  polar.}
}
\value{
  data.frame with colums:

  \item{ln.speed}{log speed estimate (in meters/time unit
  if \code{polar.coord=TRUE}.}

  \item{var.ln.speed}{log speed estimated standard error.}
}
\description{
  Calculates log speed estimates and standard errors based
  on the delta method. Typically this is used within the
  function \code{\link{crwPredict}}.
}
\author{
  Devin S. Johnson
}

