"mergeTrackStop" <- function(data, stopData, Time.name="Time",
                             interp=c('zeros','ma0'), win=2, constCol)
{
    interp <- interp[1]
    nmsStop <- names(stopData)[!names(stopData) %in% Time.name]
    if (length(nmsStop) > 1) stop("You can only merge 1 stopping variable at a time")
    stopData <- stopData[order(stopData[, Time.name], stopData[, nmsStop]), ]
    stopData <- stopData[!duplicated(stopData[, Time.name]), ]
    stopStart <- min(stopData[, Time.name][!is.na(stopData[, nmsStop])])
    stopEnd <- max(stopData[, Time.name])
    trackStart <- min(data[, Time.name])
    trackEnd <- max(data[, Time.name])
    Start <- max(stopStart, trackStart)
    End <- min(stopEnd, trackEnd)
    if (!missing(constCol)) constVal <- data[1, constCol]
    mergeData <- data.frame(seq(floor(Start), ceiling(End), 1))
    names(mergeData) <- Time.name
    stopData <- merge(stopData, mergeData, by=Time.name, all=TRUE)
    stopData$stopType <- ifelse(is.na(stopData[, nmsStop]), "p", "o")
    dtTmp <- ifelse(is.na(stopData[, nmsStop]), 0, stopData[, nmsStop])
    if (interp == "ma0") {
        stopTimeAvg <- filter(dtTmp, rep(1 / (2 * win + 1), 2 * win + 1),
                              "convolution", sides=2)
        dtTmp <- ifelse(is.na(stopData[, nmsStop]), stopTimeAvg, stopData[, nmsStop])
    }
    stopData[, paste(nmsStop, "Orig", sep="")] <- stopData[, nmsStop]
    stopData[, nmsStop] <- ifelse(is.na(dtTmp), 0, dtTmp)
    stopData$locType <- "p"
    data$locType <- "o"
    data <- merge(data, stopData, all=TRUE)
    data <- data[order(data[, Time.name]), ]
    data <- data[!duplicated(data[, Time.name]), ]
    data[, nmsStop] <- round(approx(data[, Time.name], data[, nmsStop],
                                    xout=data[, Time.name], method="constant")$y, digits=3)
    pind <- ifelse(data$stopType == "p", 1, 0)
    pind <- approx(pind, xout=1:length(pind), method="constant", rule=2)$y
    data$stopType <- ifelse(pind == 1, "p", "o")
    if (!missing(constCol)) {
        for (l in 1:length(constCol)) {
            data[, constCol[l]] <- constVal[l]
        }
    }
    out <- data[(data[, Time.name] >= trackStart & data[, Time.name] <= trackEnd), ]
    return(out)
}
