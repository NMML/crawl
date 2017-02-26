#' @title tidy-like method for crwFit object
#' @param fit \code{crwFit} object from \code{crawl::crwMLE}
#' @export
tidy_crwFit <- function(fit) {
  terms <- data.frame(term = fit$nms)
  out <- as.data.frame(round(cbind(
    fit$par, fit$se, fit$ci[, 1], fit$ci[, 2]), 3)
  )
  
  colnames(out) <- c("estimate", 
                     "std.error", 
                     "conf.low", 
                     "conf.high")
  out <- cbind(terms,out)
  
  out <- rbind(out,
               data.frame(
                 term = "logLik",
                 estimate = fit$loglik,
                 std.error = NA,
                 conf.low = NA,
                 conf.high = NA
               ),
               data.frame(
                 term = "AIC",
                 estimate = fit$aic,
                 std.error = NA,
                 conf.low = NA,
                 conf.high = NA
               ))
  out
}
