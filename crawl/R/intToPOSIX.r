`intToPOSIX` <-
function(timeVector, tz='GMT')
################################################################################
# Convert integer time to POSIX
# timeVector  = integer time in seconds
# tz          = time zone code
################################################################################
{
  Epoch = as.POSIXct(strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S",tz=tz),tz=tz)
  Epoch+timeVector
}

