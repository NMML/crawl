#' 'Flattening' a list-form crwPredict object into a data.frame
#' 
#' \dQuote{Flattens} a list form \code{\link{crwPredict}} object into a flat
#' data.frame.
#' 
#' 
#' @param predObj A crwPredict object
#' @return a \code{\link{data.frame}} version of a crwPredict list with columns
#' for the state standard errors
#' @author Devin S. Johnson
#' @seealso \code{\link{northernFurSeal}} for use example
#' @export
"as.flat" <- function(predObj)
{
    se.y <- sqrt(t(apply(predObj$V.hat.y, 3, diag)))
    se.x <- sqrt(t(apply(predObj$V.hat.x, 3, diag)))
    colnames(se.y) <- paste("se", names(predObj$alpha.hat.y), sep=".")
    colnames(se.x) <- paste("se", names(predObj$alpha.hat.x), sep=".")
    flat <- cbind(predObj$originalData, predObj$alpha.hat.y, se.y,
                  predObj$alpha.hat.x, se.x)
    if (!is.null(predObj$speed)) flat <- cbind(flat, predObj$speed)
    class(flat) <- c("crwPredict", "data.frame")
    attr(flat, "coord") <- attr(predObj, "coord")
    attr(flat, "random.drift") <- attr(predObj, "random.drift")
    attr(flat, "stop.model") <- attr(predObj, "stop.model")
    attr(flat, "polar.coord") <- attr(predObj, "polar.coord")
    attr(flat, "Time.name") <- attr(predObj, "Time.name")
    attr(flat, "flat") <- TRUE
    return(flat)
}

#' @S3method print crwFit
"print.crwFit" <- function(x, ...)
{
    fit <- x
    cat("\n\n")
    cat("Continuous-Time Correlated Random Walk fit\n\n")
    cat('Models:\n')
    cat("--------\n")
    cat("Movement   "); cat(as.character(fit$mov.model)); cat("\n")
    cat("Error   "); cat(as.character(fit$err.model)); cat("\n")
    if (fit$random.drift | !is.null(fit$stop.model)) cat("with ")
    if (fit$random.drift) cat("Random Drift")
    if (fit$random.drift & !is.null(fit$stop.model)) cat(" and ")
    if (!is.null(fit$stop.model)) cat("Movement Stops")
    cat("\n\n")
    out <- as.data.frame(round(cbind(fit$par, fit$se, fit$ci[, 1], fit$ci[, 2]), 3))
    colnames(out) <- c("Parameter Est.", "St. Err.", "95% Lower", "95% Upper")
    rownames(out) <- fit$nms
    out[!is.na(fit$fixPar), 2:4] <- "."
    print(out)
    cat("\n\n")
    cat(paste("Log Likelihood =", round(fit$loglik, 3),"\n", collapse=""))
    cat(paste("AIC =", round(fit$aic, 3),"\n", collapse=""))
    cat("\n\n\n")
}
