\name{crwMLE}
\alias{crwMLE}
\title{Fit Continuous-Time Correlated Random Walk Models to Animal Telemetry Data}
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
  covariate data. A 'SpatialPointsDataFrame' object from
  the package 'sp' will also be accepted. In which case the
  \code{polar.coord} and \code{coord} values will be taken
  from the spatial data set and ignored in the arguments.}

  \item{coord}{A 2-vector of character values giving the
  names of the "X" and "Y" coordinates in \code{data}.}

  \item{polar.coord}{logical indicating location are in
  degrees latitude and longitude.}

  \item{Time.name}{character indicating name of the
  location time column}

  \item{initial.state}{list object containg the inital
  state of the Kalman filter.}

  \item{theta}{starting values for parameter optimization.}

  \item{fixPar}{Values of parameters which are held fixed
  to the given value.}

  \item{method}{Optimization method that is passed to
  \code{\link{optim}}.}

  \item{control}{Control list which is passed to
  \code{\link{optim}}.}

  \item{constr}{Named list with elements \code{lower} and
  \code{upper} that are vectors the same length as theta
  giving the box constraints for the parameters}

  \item{prior}{A function returning the log-density
  function of the parameter prior specification}

  \item{need.hess}{A logical value which decides whether or
  not to evaluate the Hessian for parameter standard
  errors}

  \item{initialSANN}{Control list for \code{\link{optim}}
  when simulated annealing is used for obtaining start
  values. See details}

  \item{attempts}{The number of times likelihood
  optimization will be attempted}
}
\value{
  A list with the following elements:

  \item{par}{Parameter maximum likelihood estimates
  (including fixed parameters)}

  \item{estPar}{MLE without fixed parameters}

  \item{se}{Standard error of MLE}

  \item{ci}{95\% confidance intervals for parameters}

  \item{Cmat}{Parameter covariance matrix}

  \item{loglik}{Maximized log-likelihood value}

  \item{aic}{Model AIC value}

  \item{initial.state}{Intial state provided to
  \code{crwMLE} for model fitting}

  \item{coord}{Coordinate names provided for fitting}

  \item{fixPar}{Fixed parameter values provided}

  \item{convergence}{Indicator of convergence (0 =
  converged)}

  \item{message}{Meesages given by \code{optim} during
  parameter optimization}

  \item{stop.model}{Model provided for stopping variable}

  \item{random.drift}{Logical value indicating random drift
  model}

  \item{mov.model}{Model description for movement
  component}

  \item{err.model}{Model description for location error
  component}

  \item{n.par}{number of parameters}

  \item{nms}{parameter names}

  \item{n.mov}{number of movement parameters}

  \item{n.errX}{number or location error parameters for
  "longitude" error model}

  \item{n.errY}{number or location error parameters for
  "latitude" error model}

  \item{stop.mf}{covariate for stop indication in stopping
  models}

  \item{polar.coord}{Logical indicating coordinates are
  polar latitude and longitude}

  \item{init}{Initial values for parameter optimization}

  \item{data}{Original data.frame used to fit the model}

  \item{lower}{The lower parameter bounds}

  \item{upper}{The upper parameter bounds}

  \item{need.hess}{Logical value}

  \item{runTime}{Time used to fit model}
}
\description{
  The function uses the Kalman filter to estimate movement
  paramters in a state-space version of the continuous-time
  movement model. Separate models are specified for
  movement portion and the location error portion. Each
  model can depend on time indexed covariates. A
  \dQuote{haul out} model where movement is allowed to
  completely stop, as well as, a random drift model can be
  fit with this function.
}
\details{
  A full model specification involves 4 components: a
  movement model, a stopping model, 2 location error
  models, and a drift indication. The movement model
  (\code{mov.model}) specifies how the movement parameters
  should vary over time. This is a function of specified,
  time-indexed, covariates. The movement parameters (sigma
  for velocity variation and beta for velocity
  autocorrelation) are both modeled with a log link as par
  = exp(eta), where eta is the linear predictor based on
  the covariates. The \code{err.model} specification is a
  list of 2 such models, one for \dQuote{longitude} and one
  for \dQuote{latitude} (in that order) location error. If
  only one location error model is given, it is used for
  both coordinates (parameter values as well). If
  \code{drift.model} is set to \code{TRUE}, then, 2
  additional parameters are estimated for the drift
  process, a drift variance and a beta multiplier. If
  \code{polar.coord=TRUE} then the ad-hoc logitude
  correction factor described by Johnson et al. (2008)
  (Ecology 89:1208-1215) is used to adjust the variance
  scale for the longitude mdoel.

  The \code{inital.state} is a list with the following
  elemets (with the exact names):

  \code{a1.y} A vector with initial state values for the
  \dQuote{latitude} coordinate. It have 2 elemets (location
  at time 1, velocity at time 1) for non-drift models and 3
  elemets for drift models (location at time 1, velocity at
  time 1, drift velocity at time 1) for driftmodels,

  \code{P1.y} Covarince matrix for the state at time 1
  (measure of uncertainty for your inital state)
  \code{a1.y},

  \code{a1.x} Same as \code{a1.y}, but in the
  \dQuote{longitude} coordinate,

  \code{P1.x} Same as \code{P1.y}, but in the
  \dQuote{longitude} coordinate.

  \code{theta} and \code{fixPar} are vectors with the
  appropriate number or parameters. \code{theta} contains
  only those paraemters which are to be estimated, while
  \code{fixPar} contains all parameter values with
  \code{NA} for parameters which are to be estimated.

  The data set specified by \code{data} must contain a
  numeric or POSIXct column which is used as the time index
  for analysis. The column name is specified by the
  \code{Time.name} argument. If a POSIXct column is used it
  is internally converted to a numeric vector with units of
  hours. If this is not appropriate for your data, it is
  better to convert it yourself prior to analysis with
  crawl. Also, for stopping models, the stopping covariate
  must be between 0 and 1 inclusive, with 1 representing
  complete stop of the animal (no true movement, however,
  location error can still occur) and 0 represent
  unhindered movement. The coordinate location should have
  \code{NA} where no location is recorded, but there is a
  change in the movment covariates.

  The CTCRW models can be difficult to provide good initial
  values for optimization. If \code{initialSANN} is
  specified then simulated annealing is used first to
  obtain starting values for the specified optimaization
  method. If simulated annealing is used first, then the
  returned \code{init} list of the crwFit object will be a
  list with the results of the simulated annealing
  optimization.
}
\examples{
data(harborSeal)
head(harborSeal)
## Calculate Log multipliers for Argos error
argosClasses <- c("3", "2", "1", "0", "A", "B")
ArgosMultFactors <- data.frame(Argos_loc_class=argosClasses,
errX=log(c(1, 1.5, 4, 14, 5.21, 20.78)),
errY=log(c(1, 1.5, 4, 14, 11.08, 31.03)))
hsNew <- merge(harborSeal, ArgosMultFactors, by=c("Argos_loc_class"), all=TRUE)
hsNew <- hsNew[order(hsNew$Time), ]
head(hsNew)

## Initial state values
initial.dry <- list(
a1.x=c(harborSeal$longitude[1],0),
a1.y=c(harborSeal$latitude[1],0),
P1.x=diag(c(1,1)),
P1.y=diag(c(1,1))
)

##Fit model as given in Johnson et al. (2008) Ecology 89:1208-1215
## Start values for theta come from the estimates in Johnson et al. (2008)

fit1 <- crwMLE(
mov.model=~1, err.model=list(x=~errX, y=~errY), stop.model=~DryTime,
data=hsNew, coord=c("longitude","latitude"), polar.coord=TRUE, Time.name="Time",
initial.state=initial.dry, fixPar=c(NA, 1, NA, 1, NA, NA, NA), theta=c(-6,-7,-4,-0.5,-1),
control=list(maxit=2000, trace=1, REPORT=1),
)

fit1
str(fit1)

##Use simulated annealing to obtain start values and place constraints on the parameters

fit2 <- crwMLE(
mov.model=~1, err.model=list(x=~errX, y=~errY), stop.model=~DryTime,
data=hsNew, coord=c("longitude","latitude"), Time.name="Time", polar.coord=TRUE,
initial.state=initial.dry, fixPar=c(NA, 1, NA, 1, NA, NA, NA),
control=list(maxit=2000, trace=1, REPORT=1),
constr=list(lower=c(-6, -6, -Inf, -Inf, -Inf), upper=Inf),
initialSANN=list(maxit=100, temp=5, tmax=5, trace=1, REPORT=2)
)

fit2

##See simulated annealing start values
fit2$init$par
}
\author{
  Devin S. Johnson
}
\seealso{
  \code{\link{northernFurSeal}} for additional examples.
}

