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
