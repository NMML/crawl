#' Detect appropriate time scale for movement analysis
#' 
#' This function examines the time vector and evaluates the median time
#' interval. With this, we determine what the best time scale for the 
#' movement model is likely to be.
#'
#' @param time_vector a vector of class POSIXct
#'
#' @return character of either "seconds","minutes","hours","days","weeks"
#' @export
#'

detect_timescale <- function(time_vector) {
  intervals <- difftime(time_vector[-length(time_vector)], 
                        time_vector[-1], units = "secs")
  median_int <- abs(median(intervals))
  if (median_int < 31) {
    return("seconds")
  }
  if (median_int < 1801) {
    return("minutes")
  }
  if (median_int < 3600*12) {
    return("hours")
  }
  if (median_int < 3600*24*3.5) {
    return("days")
  }
  return("weeks")
}
