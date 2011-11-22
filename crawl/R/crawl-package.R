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
#' 1.3-2\cr Date: \tab 2011-11-21\cr License: \tab Unlimited \cr LazyLoad: \tab
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
#' @examples
#' 
#' 
#' data(northernFurSeal)
#' 
#' argosClasses <- c("3", "2", "1", "0", "A", "B")
#' ArgosMultFactors <- data.frame(Argos_loc_class=argosClasses,
#'                                errX=log(c(1, 1.5, 4, 14, 5.21, 20.78)),
#'                                errY=log(c(1, 1.5, 4, 14, 11.08, 31.03)))
#' nfsNew <- merge(northernFurSeal, ArgosMultFactors,
#'                 by=c("Argos_loc_class"), all.x=TRUE)
#' nfsNew <- nfsNew[order(nfsNew$Time), ]
#' 
#' # State starting values
#' initial.drift <- list(a1.x=c(189.686, 0, 0), a1.y=c(57.145, 0, 0),
#'                       P1.x=diag(c(0, 0.001, 0.001)),
#'                       P1.y=diag(c(0, 0.001, 0.001)))
#' 
#' ##Fit random drift model
#' fit <- crwMLE(mov.model=~1, err.model=list(x=~errX, y=~errY), drift.model=TRUE,
#'               data=nfsNew, coord=c("longitude", "latitude"), polar.coord=TRUE,
#'               Time.name="Time", initial.state=initial.drift, 
#'               fixPar=c(NA, 1, NA, 1, NA, NA, NA, NA), 
#'               control=list(maxit=2000,trace=1, REPORT=10),
#'               initialSANN=list(maxit=300, trace=1, REPORT=1)
#'               )
#' 
#' ##Make hourly location predictions
#' predTime <- seq(ceiling(min(nfsNew$Time)), floor(max(nfsNew$Time)), 1)
#' predObj <- crwPredict(object.crwFit=fit, predTime, speedEst=TRUE, flat=TRUE)
#' head(predObj)
#' crwPredictPlot(predObj)
#' 
#' ##Create simulation object with 100 parameter draws
#' simObj <- crwSimulator(fit, predTime, parIS=100, df=20, scale=18/20)
#' 
#' ## Examine IS weight distribution
#' w <- simObj$thetaSampList[[1]][,1]
#' dev.new()
#' hist(w*100, main='Importance Sampling Weights', sub='More weights near 1 is desirable')
#' 
#' ##Approximate number of independent samples
#' round(100/(1+(sd(w)/mean(w))^2))
#' 
#' dev.new(bg=gray(0.75))
#' jet.colors <-
#'   colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#'                      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#' crwPredictPlot(predObj, 'map')
#' 
#' ## Sample 20 tracks from posterior predictive distribution
#' iter <- 20
#' cols <- jet.colors(iter)
#' for(i in 1:iter){
#'    samp <- crwPostIS(simObj)
#'    lines(samp$alpha.sim.x[,'mu'], samp$alpha.sim.y[,'mu'],col=cols[i])
#' }
#' 
#' 
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



.onLoad <- function(library, pkgname)
{
  ## Return a list, each element of which is a vector
  ## the first element of the vector is the stuff before the colon in info[[1]]
  ## the second element is the stuff after the colon (can get > 2 elements some
  ## times but ignore)
  info <- strsplit(library(help=pkgname, character.only=TRUE)$info[[1]], "\\:[ ]+")
  ## Go through the list, pulling out the Package, Version and Built strings
  l <- length(info)
  package <- version <- built <- ""
  for (i in 1:l) {
    if(info[[i]][1] == "Package") package <- info[[i]][2]
    if(info[[i]][1] == "Version") version <- info[[i]][2]
    if(info[[i]][1] == "Built") built <- info[[i]][2]
  }
  ## Print these out
  packageStartupMessage(paste("This is", package, version, "\nBuilt:", built, "\n"))
  ## uncomment for fortran/c code
  ## library.dynam("filenameForDll", pkgname)
  #library.dynam("crawl", pkgname)
}

.onUnload <- function(libpath)
{
  #library.dynam.unload("crawl", libpath)
  cat("\nBye-Bye from crawl\n\n")
  return(invisible())
}

