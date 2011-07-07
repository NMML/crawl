\name{fillCols}

\alias{fillCols}

% -------------------------------------------------------------------------

\title{Fill missing values in data set (or matrix) columns for which there is a
       single unique value}


\description{

  Looks for columns in a data set that have a single unique non-missing value
  and fills in all \code{NA} with that value

}


\usage{

fillCols(data)

}


\arguments{

  \item{data}{data.frame}

}


\value{ data.frame }

% -------------------------------------------------------------------------

\examples{
#library(crawl)
data1 <- data.frame(constVals=rep(c(1,NA),5), vals=1:10)
data1[5,2] <- NA
data1
data2 <- fillCols(data1)
data2

mat1 <- matrix(c(rep(c(1,NA),5), 1:10), ncol=2)
mat1[5,2] <- NA
mat1
mat2 <- fillCols(mat1)
mat2


}


\author{Devin S. Johnson}
