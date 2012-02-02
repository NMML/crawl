\name{as.flat}
\alias{as.flat}
\title{'Flattening' a list-form crwPredict object into a data.frame}
\usage{
  as.flat(predObj)
}
\arguments{
  \item{predObj}{A crwPredict object}
}
\value{
  a \code{\link{data.frame}} version of a crwPredict list
  with columns for the state standard errors
}
\description{
  \dQuote{Flattens} a list form \code{\link{crwPredict}}
  object into a flat data.frame.
}
\author{
  Devin S. Johnson
}
\seealso{
  \code{\link{northernFurSeal}} for use example
}

