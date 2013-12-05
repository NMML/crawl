\name{displayPar}
\alias{displayPar}
\title{Display the order of parameters along with fixed values and starting values}
\usage{
  displayPar(mov.model = ~1, err.model = NULL,
    stop.model = NULL, drift.model = FALSE, data, theta,
    fixPar)
}
\arguments{
  \item{mov.model}{formula object specifying the time
  indexed covariates for movement parameters.}

  \item{err.model}{A 2-element list of formula objects
  specifying the time indexed covariates for location error
  parameters.}

  \item{stop.model}{formula object giving the covariate for
  the stopping portion of the model.}

  \item{drift.model}{logical indicating whether or not to
  include a random drift component.}

  \item{data}{data.frame object containg telemetry and
  covariate data. A \code{\link{SpatialPointsDataFrame}}
  object from the package 'sp' will also be accepted.}

  \item{theta}{starting values for parameter optimization.}

  \item{fixPar}{Values of parameters which are held fixed
  to the given value.}
}
\value{
  A data frame with the following columns

  \item{ParNames}{The names of the parameters specified by
  the arguments.}

  \item{fixPar}{The values specified by the \code{fixPar}
  argument for fixed values of the parameters. In model
  fitting, these values will remain fixed and will not be
  estimated.}

  \item{thetaIndex}{This column provides the index of each
  element of the theta argument and to which parameter it
  corresponds.}

  \item{thetaStart}{If a value is given for the
  \code{theta} argument it will be placed in this column
  and its elements will correspond to the \code{thetaIdx}
  column.}
}
\description{
  This function takes the model spesification arguments to
  the \code{\link{crwMLE}} function and displays a table
  with the parameter names in the order that \code{crwMLE}
  will use during model fitting. This is useful for
  specifying values for the \code{fixPar} or \code{theta}
  (starting values for free parameters) arguments.
}
\author{
  Devin S. Johnson
}
\seealso{
  \code{demo(northernFurSealDemo)} for example.
}

