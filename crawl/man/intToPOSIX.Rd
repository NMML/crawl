\name{intToPosix}

\alias{intToPOSIX}

% -------------------------------------------------------------------------

\title{Reverse as.numeric command that is performed on a vector of type POSIXct}


\description{

  Takes integer value produced by \code{as.numeric(x)}, where \code{x} is a POSIXct
  vector and returns it to a POSIXct vector

}


\usage{

intToPOSIX(timeVector, tz = "GMT")

}


\arguments{

  \item{timeVector}{A vector of integers produced by as.numeric applied to a PSIXct vector}
  
  \item{tz}{Time zone of the vector (see \code{\link{as.POSIXct}}).}

}


\value{POSIXct vector}

\note{
   There is no check that as.numeric applied to a POSIX vector produced \code{timeVector}.
   So, caution is required in using this function. It was included simply
   because I have found it useful
}


% -------------------------------------------------------------------------

\examples{
#library(crawl)
timeVector <- as.numeric(Sys.time())
timeVector
intToPOSIX(timeVector, tz="")

}


\author{Devin S. Johnson}
