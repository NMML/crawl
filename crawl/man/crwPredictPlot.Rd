\name{crwPredictPlot}

\alias{crwPredictPlot}

% -------------------------------------------------------------------------

\title{Plot CRW predicted object}


\description{

  Creates 2 types of plots of a crwPredict object: a plot of both
  coordinate axes with prediction intervals and a plot of just observed
  locations and predicted locations.

}


\usage{
crwPredictPlot(object, plotType="ll")

}


\arguments{

  \item{object}{\code{crwPredict} object.}

  \item{plotType}{type of plot has to be one of the following:
    \dQuote{map} or \dQuote{ll} (default).}

}


\value{ A plot. }


\seealso{ See \code{\link{northernFurSeal}} for additional examples. }

% -------------------------------------------------------------------------

\author{Devin S. Johnson and Sebastian Luque}
