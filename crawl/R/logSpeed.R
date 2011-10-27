#' calculate estimate of log speed from a crwPredict object
#' 

#' 
#' Calculates log speed estimates and standard errors based on the delta
#' method. Typically this is used within the function \code{\link{crwPredict}}.
#' 
#' 
#' @param predx matrix of longitude state estimates.
#' @param predy matrix of latitude state estimates.
#' @param varx an array of covariance matrices for predx.
#' @param vary ana array of covariance matrices for predy.
#' @param polar.coord Logical. TRUE if coordinates are polar.
#' @return
#' 
#' data.frame with colums:
#' 
#' \item{ln.speed}{log speed estimate (in meters/time unit if
#' \code{polar.coord=TRUE}.}
#' 
#' \item{var.ln.speed}{log speed estimated standard error.}
#' @author Devin S. Johnson
#' @export
"logSpeed" <- function(predx, predy, varx, vary, polar.coord)
{
    rd <- ncol(predx) == 3
    pc <- polar.coord
    if (pc) {
        lnspd <- expression(0.5 * log((111325 * (ay2 + ay3)) ^ 2 +
            (111325 * (ax2 + ax3)) ^ 2 * cos(ay1 * pi / 180) ^ 2))
    } else {
        lnspd <- expression(0.5 * log((ay2 + ay3) ^ 2 + (ax2 + ax3) ^ 2))
    }
    d.lnspd.x <- deriv(lnspd, c("ax1", "ax2", "ax3"),
                       function.arg=c("ax1", "ax2", "ax3", "ay1", "ay2", "ay3"))
    d.lnspd.y <- deriv(lnspd,c("ay1", "ay2", "ay3"),
                       function.arg=c("ax1", "ax2", "ax3", "ay1", "ay2", "ay3"))
    if (rd) {
        ax3 <- predx[, 3]
        ay3 <- predy[, 3]
    } else {
        ax3 <- 0
        ay3 <- 0
    }
    der.a.x <- attr(d.lnspd.x(predx[, 1], predx[, 2], ax3,
                              predy[, 1], predy[, 2], ay3), "gradient")
    der.a.y <- attr(d.lnspd.y(predx[, 1], predx[, 2], ax3,
                              predy[, 1], predy[, 2], ay3), "gradient")
    ln.speed <- d.lnspd.x(predx[, 1], predx[, 2], ax3, predy[, 1], predy[, 2], ay3)
    var.lnspd <- sapply(1:nrow(predx), function(i) {
        t(der.a.x[i,1:(2 + rd)]) %*% varx[, , i] %*% der.a.x[i, 1:(2 + rd)] +
                       t(der.a.y[i, 1:(2 + rd)]) %*% vary[, , i] %*%
                           der.a.y[i, 1:(2 + rd)]
    })
    return(data.frame(ln.speed=as.double(ln.speed), var.ln.speed=as.double(var.lnspd)))
}
