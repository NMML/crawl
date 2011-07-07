\name{harborSeal}
\docType{data}

\alias{harborSeal}

% -------------------------------------------------------------------------

\title{Harbor seal telemetry data}


\description{

  Harbor seal relocation data set used in Johnson et al. (2008)

}


\usage{

data(harborSeal)

}


\format{

  A data frame with 7059 observations on the following 5 variables.

  \describe{
    \item{\code{Time}}{a numeric vector.}

    \item{\code{latitude}}{a numeric vector.}

    \item{\code{longitude}}{a numeric vector.}

    \item{\code{DryTime}}{a numeric vector.}

    \item{\code{Argos_loc_class}}{a factor with levels \code{0} \code{1}
      \code{2} \code{3} \code{A} \code{B}}.}

}


\source{
  Polar Ecosystems Program
  National Marine Mammal Laboratory
  Alaska Fisheries Science Center
  National Marine Fisheries Service, NOAA
  7600 Sand Point Way, NE
  Seattle, WA 98115
}


\references{ Johnson, D., J. London, M. -A. Lea, and J. Durban (2008)
Continuous-time random walk model for animal telemetry data. Ecology
89:1208-1215. }

% -------------------------------------------------------------------------

\examples{

data(harborSeal)
head(harborSeal)

}


\author{Devin S. Johnson}


\keyword{datasets}
