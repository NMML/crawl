#' calculate estimate of log speed from a crwPredict object
#' 

#' 
#' Calculates log speed estimates and standard errors based on the delta
#' method. Typically this is used within the function \code{\link{crwPredict}}.
#' 
#' 
#' @param alpha.hat State estimates from \code{\link{crwPredict}}
#' @param V.hat An array of covariance matrices for alpha.hat.
#' @return
#' 
#' data.frame with colums:
#' 
#' \item{ln.speed}{log speed estimate (in meters/time unit)}
#' \item{var.ln.speed}{log speed estimated standard error.}
#' @author Devin S. Johnson
#' @export
"logSpeed" <- function(alpha.hat, V.hat)
{
    rd <- ncol(predx) == 3
    lnspd <- expression(0.5 * log((ay2 + ay3) ^ 2 + (ax2 + ax3) ^ 2))
    d.lnspd.x <- deriv(lnspd, c("ax1", "ax2", "ax3"),
                       function.arg=c("ax1", "ax2", "ax3", "ay1", "ay2", "ay3"))
    d.lnspd.y <- deriv(lnspd,c("ay1", "ay2", "ay3"),
                       function.arg=c("ax1", "ax2", "ax3", "ay1", "ay2", "ay3"))
    if (rd) {
        ax3 <- alpha.hat[, 3]
        ay3 <- alpha.hat[, 6]
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
