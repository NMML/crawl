#' Generic subset/bracket method for crwIS classes
#' 
#' @param ... other arguments
#' @export

"[.crwIS" <- function(x, i, j, ... , drop = TRUE)
          {
            x$alpha.sim <- x$alpha.sim[i, ]
            x$locType <- x$locType[i]
            x$TimeNum <- x$TimeNum[i]
            x$Time <- x$Time[i]
            return(x)
          }
