#' Coerce crawl objects (crwIS and crwPredict) to tibbles
#'
#' @author Josh M. London
#' @param crw_object an object of class \code{"crwIS"} or \code{"crwPredict"}
#' @param ... Additional arguments that are ignored
#' @export

crw_as_tibble <- function(crw_object, ...) {
  UseMethod("crw_as_tibble",crw_object)
}

#' @describeIn crw_as_tibble coerce crwIS object to tibble
#' @export
crw_as_tibble.crwIS <- function(crw_object, ...) {
  tn = attr(crw_object, "Time.name")
  out = data.frame(TimeNum=crw_object$TimeNum, locType =crw_object$locType, crw_object$alpha.sim)
  out[,tn] = crw_object[[tn]]
  out=tibble::as_tibble(out) %>% dplyr::arrange(.data$TimeNum)
  out
}

#' @describeIn crw_as_tibble coerce crwPredict object to tibble
#' @export
crw_as_tibble.crwPredict <- function(crw_object, ...) {
  if(inherits(crw_object,"list")){
    crw_object = fillCols(crawl::flatten(crw_object))
  }
  tibble::as_tibble(crw_object) 
}

#' @describeIn crw_as_tibble 
#' @export
crw_as_tibble.tbl <- function(crw_object, ...) {
  crw_object
}
