library(matrixStats)
data("pt01EcoG")
timeWindow <- c(-30, 20)
epoch <- Epoch(pt01EcoG)
fs=1000
sozIndex <- attr(pt01EcoG, "sozIndex")
timeNum <- ncol(epoch)
windowParams<-c(0.25,0.1)


# Set spectrogram parameters
# define frequency bands
deltaBand<-c(0.5,4)
thetaBand<-c(4,8)
alphaBand<-c(8,13)
betaBand<-c(13,30) 
gammaBand<-c(30,90) 
highGammaBand<-c(90,150)

rangeBand<-betaBand

nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1
data   <- vector(mode="numeric", length=timeNum)
data[1:timeNum]<-epoch$data[sozIndex[1],1:timeNum]
# Compute the multitaper spectrogram
results = multitaperSpectrogramR(data=data, fs=fs, windowParams = windowParams, frequencyRange=rangeBand)
## Visualize a subject of electrodes
#display <- c(sozIndex, 77:80)
stimes = results[[2]]

timesOnset<-results[[2]]+timeWindow[1]

endBaseIndex<-which.min(abs(timesOnset + 20))

  

elecNum <- nrow(epoch)
nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1

ermaster=matrix(0,elecNum,nwt)

# comute one user defined band

for(ie in 1:elecNum){
  
  data[1:timeNum]<-epoch$data[ie,1:timeNum]
  # Compute the multitaper spectrogram
  results = multitaperSpectrogramR(data, fs)
  spec=results[[1]]
  sfreq=results[[3]]
  
  
  freqStart=deltaBand[1]
  freqEnd=deltaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectDelta = results[[1]][freqStartIndex:freqEndIndex,]
  powDelta=colSums(spectDelta)
  mBaseDelta<-
  
  freqStart=thetaBand[1]
  freqEnd=thetaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  specttheta = results[[1]][freqStartIndex:freqEndIndex,]
  powtheta=colSums(specttheta)
  
  freqStart=alphaBand[1]
  freqEnd=alphaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectalpha = results[[1]][freqStartIndex:freqEndIndex,]
  powalpha=colSums(spectalpha)
  
  freqStart=betaBand[1]
  freqEnd=betaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectbeta = results[[1]][freqStartIndex:freqEndIndex,]
  powbeta=colSums(spectbeta)
  
  freqStart=gammaBand[1]
  freqEnd=gammaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectgamma = results[[1]][freqStartIndex:freqEndIndex,]
  powgamma=colSums(spectgamma)

  
  freqStart=highGammaBand[1]
  freqEnd=highGammaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  specthighGamma = results[[1]][freqStartIndex:freqEndIndex,]
  powhighGamma=colSums(specthighGamma)
  
}
