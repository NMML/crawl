#' Generic subset/bracket method for crwIS classes
#' @param x crwIS object
#' @param i elements to extract or replace. These are numeric or character or, 
#' empty or logical. Numeric values are coerced to integer as if by \code{as.integer}
#' @param ... other arguments
#' @param drop logical. If TRUE the result is coerced to the lowest possible 
#' dimension. 
#' @export

"[.crwIS" <- function(x, i, ..., drop = TRUE) {
  x$alpha.sim <- x$alpha.sim[i,drop=drop]
  x$locType <- x$locType[i,drop=drop]
  x$TimeNum <- x$TimeNum[i,drop=drop]
  x[[attr(x, "Time.name")]] <- x[[attr(x, "Time.name")]][i,drop=drop]
  return(x)
}
