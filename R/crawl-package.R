#' @title Fit Continuous-Time Correlated Random Walk Models to Animal Movement Data
#' 
#' @description The Correlated RAndom Walk Library (I know it is not an R library,
#' but, "crawp" did not sound as good) of R functions was designed for fitting
#' continuous-time correlated random walk (CTCRW) models with time indexed
#' covariates. The model is fit using the Kalman-Filter on a state space
#' version of the continuous-time staochistic movement process.
#' 
#' \tabular{ll}{ 
#' Package: \tab crawl\cr 
#' Type: \tab Package\cr 
#' Version: \tab 2.2.1\cr 
#' Date: \tab September 12, 2018\cr 
#' License: \tab CC0 \cr 
#' LazyLoad: \tab yes\cr 
#' }
#' 
#' @note This software package is developed and maintained by scientists at the NOAA Fisheries Alaska 
#' Fisheries Science Center and should be considered a fundamental research communication. 
#' The reccomendations and conclusions presented here are those of 
#' the authors and this software should not be construed as official communication by NMFS, NOAA, 
#' or the U.S. Dept. of Commerce. In addition, reference to trade names does not imply endorsement by the 
#' National Marine Fisheries Service, NOAA. While the best efforts have been made to insure the 
#' highest quality, tools such as this are under constant development and are subject to change.
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
#' @import dplyr
#' @importFrom Rcpp evalCpp
#' @importFrom graphics layout
#' @importFrom methods as slot 
#' @importFrom stats approx model.frame model.matrix 
#'             na.pass optim pchisq pexp pnorm qnorm 
#'             rchisq runif sd setNames
#' @useDynLib crawl, .registration = TRUE

NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Northern fur seal pup relocation data set used in Johnson et al. (2008)
#' 
#' 
#' @name northernFurSeal
#' @docType data
#' @format A data frame with 795 observations on the following 4 variables:
#' 
#' \describe{ \item{GMT}{A POSIX time vector}
#' 
#' \item{loc_class}{a factor with levels \code{3} \code{2}
#' \code{1} \code{0} \code{A}.}
#' 
#' \item{lat}{a numeric vector. Latitude for the locations}
#' 
#' \item{long}{a numeric vector. Longitude for the locations}
#' 
#' }
#' @references Johnson, D., J. London, M. -A. Lea, and J. Durban (2008) Continuous-time
#' random walk model for animal telemetry data. Ecology 89:1208-1215.
#' @source Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA 7600 Sand
#' Point Way NE Seattle, WA 98115
#' @keywords datasets
NULL

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
#' @source Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA 7600 Sand
#' Point Way NE Seattle, WA 98115
#' @keywords datasets
NULL

#' Bearded Seal Location Data
#' 
#' 
#' @name beardedSeals
#' @docType data
#' @format A data frame with 27,548 observations on 3 bearded seals in Alaska:
#' 
#' \describe{
#' \item{deployid}{Unique animal ID}
#' \item{ptt}{Hardware ID}
#' \item{instr}{Hardware type}
#' \item{date_time}{Time of location}
#' \item{type}{Location type}
#' \item{quality}{Argos location quality}
#' \item{latitude}{Observed latitude}
#' \item{longitude}{Observed longitude}
#' \item{error_radius}{Argos error radius}
#' \item{error_semimajor_axis}{Argos error ellipse major axis length}
#' \item{error_semiminor_axis}{Argos error ellipse minor axis length}
#' \item{error_ellipse_orientation}{Argos error ellispse degree orientation}
#' }
#' @source Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA 7600 Sand
#' Point Way NE Seattle, WA 98115
#' @keywords datasets
NULL




.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
    paste(paste(package, version, paste("(",date, ")", sep=""), "\n"), 
          "Demos and documentation can be found at our new GitHub repository:\n",
          "https://dsjohnson.github.io/crawl_examples/")
  )
}

# .onUnload <- function(libpath)
# {
#   #library.dynam.unload("crawl", libpath)
#   cat("\nBye-Bye from crawl\n\n")
#   return(invisible())
# }

