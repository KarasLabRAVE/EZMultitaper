standardizeIEEG <- function(data) {
  scaling <- 10^floor(log10(max(data)))
  plotData <- data / scaling
}

meanPowBaselineBand <- function(Band,spect, index){

    
  
}
