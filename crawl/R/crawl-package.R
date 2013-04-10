#' Fit Continuous-Time Correlated Random Walk models to animal movement data
#' 
#' The (C)orrelated (RA)ndom (W)alk (L)ibrary (I know it is not an R library,
#' but, "crawp" did not sound as good) of R functions was designed for fitting
#' continuous-time correlated random walk (CTCRW) models with time indexed
#' covariates. The model is fit using the Kalman-Filter on a state space
#' version of the continuous-time staochistic movement process.
#' 
#' 
#' \tabular{ll}{ Package: \tab crawl\cr Type: \tab Package\cr Version: \tab
#' 1.4\cr Date: \tab February 1, 2013\cr License: \tab Unlimited \cr LazyLoad: \tab
#' yes\cr }
#' 
#' @name crawl-package
#' @aliases crawl-package crawl
#' @docType package
#' @author Devin S. Johnson
#' 
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#' @references Johnson, D., J. London, M. -A. Lea, and J. Durban (2008)
#' Continuous-time correlated random walk model for animal telemetry data.
#' Ecology 89(5) 1208-1215.
NULL

#' Northern fur seal pup relocation data set used in Johnson et al. (2008)
#' 
#' 
#' @name northernFurSeal
#' @docType data
#' @format A data frame with 795 observations on the following 4 variables:
#' 
#' \describe{ \item{Time}{a numeric vector.}
#' 
#' \item{Argos_loc_class}{a factor with levels \code{0} \code{1}
#' \code{2} \code{3} \code{A}.}
#' 
#' \item{latitude}{a numeric vector.}
#' 
#' \item{longitude}{a numeric vector.}
#' 
#' }
#' @references
#' 
#' Johnson, D., J. London, M. -A. Lea, and J. Durban (2008) Continuous-time
#' random walk model for animal telemetry data. Ecology 89:1208-1215.
#' @source Alska Ecosystems Program National Marine Mammal Laboratory Alaska
#' Fisheries Science Center National Marine Fisheries Service, NOAA 7600 Sand
#' Point Way NE Seattle, WA 98115
#' @keywords datasets
#' @useDynLib crawl

#' Harbor seal relocation data set used in Johnson et al. (2008)
#' 
#' 
#' @name harborSeal
#' @docType data
#' @format
#' 
#' A data frame with 7059 observations on the following 5 variables.
#' 
#' \describe{ \item{Time}{a numeric vector.}
#' 
#' \item{latitude}{a numeric vector.}
#' 
#' \item{longitude}{a numeric vector.}
#' 
#' \item{DryTime}{a numeric vector.}
#' 
#' \item{Argos_loc_class}{a factor with levels \code{0} \code{1}
#' \code{2} \code{3} \code{A} \code{B}}.}
#' @author Devin S. Johnson
#' @references Johnson, D., J. London, M. -A. Lea, and J. Durban (2008)
#' Continuous-time random walk model for animal telemetry data. Ecology
#' 89:1208-1215.
#' @source Polar Ecosystems Program National Marine Mammal Laboratory Alaska
#' Fisheries Science Center National Marine Fisheries Service, NOAA 7600 Sand
#' Point Way, NE Seattle, WA 98115
#' @keywords datasets
#' @examples
#' 
#' 
#' data(harborSeal)
#' head(harborSeal)
NULL



.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
    paste(paste(package, version, paste("(",date, ")", sep=""), "\n"), 
          "Type 'demo(package='crawl')' to see a list of demos for this package.\n",
          "The raw code for the demos can be found by typing 'system.file('demo', package='crawl')'")
                        )

}

# .onUnload <- function(libpath)
# {
#   #library.dynam.unload("crawl", libpath)
#   cat("\nBye-Bye from crawl\n\n")
#   return(invisible())
# }

