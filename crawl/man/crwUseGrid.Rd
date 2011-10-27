\name{crwUseGrid}
\alias{crwUseGrid}
\title{Compute a spatial use grid from a crawl prediction}
\usage{
  crwUseGrid(object, grid, rm.zeros = FALSE, subset = TRUE)
}
\arguments{
  \item{object}{A crwPredict object output from
  \code{\link{crwPredict}} or crwIS object from
  \code{\link{crwPostIS}}.  The object MUST have been
  created using a set of prediction times (i.e.,
  is.null(predTime)=FALSE)}

  \item{grid}{A \code{GridTopology} object from the 'sp'
  package}

  \item{rm.zeros}{Logical. If set to true, the zeros in the
  use grid are set to NA for better plots. Use grids can
  have an overwhelming number of zeros.}

  \item{subset}{An indicator of which times should be used
  for calculation of the use grid. Can be a logical vector
  or a vector of integers, such as from a call to
  \code{which}}
}
\value{
  A \code{SpatialGridDataFrame} with data column 'use'
  which gives the count of predicted (see the 'sp' package)
  locations within the grid cell
}
\description{
  This function take a crwPredict or crwIS object and a
  spatial GridTopology object from the 'sp' package and
  outputs the number of predicted locations in each grid
  cell
}
\author{
  Devin S. Johnson <devin.johnson@noaa.gov>
}

