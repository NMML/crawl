"expandPred" <- function(x, Time='Time', predTime, time.col=FALSE)
{
   if(is.character(Time)) {
      Time.name <- Time
      if(!Time.name%in%colnames(x)) stop(paste(Time.name, 'is not in x. PLease specify correct time indicator'))
   }
   else if(is.numeric(Time) & length(Time)==nrow(as.matrix(x))) {
              x <- cbind(Time=Time, x)
              Time.name <- 'Time'
           }
   else stop("Value given for 'Time' is not a recognized format. See crawl documentation")
   predData <- data.frame(predTime)
   colnames(predData) <- Time.name
   newx <- merge(as.data.frame(x), predData,
                 by=c(Time.name), all=TRUE)
   for(i in 1:ncol(newx)) {
      vec <- newx[,i]
      newx[,i] <- vec[!is.na(vec)][cumsum(!is.na(vec))]
   }
   
   if(time.col) return(newx[!duplicated(newx[,Time.name]),])
   else return(newx[!duplicated(newx[,Time.name]),!colnames(newx)%in%Time.name])
}
