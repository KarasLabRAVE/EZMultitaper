standardizeIEEG <- function(data) {
  scaling <- 10^floor(log10(max(data)))
  plotData <- data / scaling
}

#' Title
#'
#' @param epoch 
#' @param fs 
#' @param windowParams 
#' @param rangeBand 
#' @param powTimeWindow 
#' @param baseTimeWindow 
#'
#' @return
#' @export
#'
#' @examples
meanPowBaselineBand <- function(epoch, fs, windowParams, rangeBand, powTimeWindow, baseTimeWindow){

  timeNum <- ncol(epoch)
  nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1
  data   <- vector(mode="numeric", length=timeNum)
  data[1:timeNum]<-epoch$data[1,1:timeNum]
  # Compute the multitaper spectrogram
  results = multitaperSpectrogramR(data=data, fs=fs, windowParams = windowParams, frequencyRange=rangeBand)
  stimes = results[[2]]
  
  timesOnset<-results[[2]]+timeWindow[1]

  startBaseIndex<-which.min(abs(timesOnset - baseTimeWindow[1]))
  endBaseIndex<-which.min(abs(timesOnset - baseTimeWindow[2]))

  elecNum <- nrow(epoch)
  
  startEpochIndex<-which.min(abs(timesOnset - powTimeWindow[1]))
  endEpochIndex<-which.min(abs(timesOnset - powTimeWindow[2]))

  timesOnset<-timesOnset[startEpochIndex:endEpochIndex]
  powBandMaster=matrix(0,elecNum,length(timesOnset))
  
  
  for(ie in 1:elecNum){
    
    data[1:timeNum]<-epoch$data[ie,1:timeNum]
    # Compute the multitaper spectrogram
    results <-multitaperSpectrogramR(data=data, fs=fs, windowParams = windowParams, frequencyRange=rangeBand)
    spect<-results[[1]]
    
    powBandOneElec<-colSums(spect)
    
    mBasePow<-mean( powBandOneElec[startBaseIndex:endBaseIndex])
    powBandOneElec<- powBandOneElec-mBasePow
    
    powBandMaster[ie,1:length(timesOnset)]<-powBandOneElec[startEpochIndex:endEpochIndex]
    
  }
  
  rownames(powBandMaster)<-rownames(epoch)
  colnames(powBandMaster)<-timesOnset
  
  powBandMaster
  
  
}
