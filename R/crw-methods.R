#' Generic methods for crwIS and crwPredict classes
#' 
#' @param drop_obs boolean whether to drop the "o" locType records
#' @param ... other arguments
#' @export

# "[.crwIS" <- function(x, i, j, drop_obs = FALSE, ... , drop = TRUE)
#           {
#             alpha.sim <- x$alpha.sim
#             locType <- x$locType
#             Time <- x$Time
#             
#             if (drop_obs) {
#               p_idx <- locType == "p"
#               alpha.sim <- alpha.sim[p_idx, ]
#               locType <- locType[p_idx]
#               Time <- Time[p_idx]
#             }
#             
#             alpha.sim <- alpha.sim[i, ]
#             locType <- locType[i]
#             Time <- Time[i]
#             
#             x$alpha.sim <- alpha.sim
#             x$locType <- locType
#             x$Time <- Time
#             
#             return(x)
#           }
