library(matrixStats)
data("pt01EcoG")
timeWindow <- c(-30, 20)
epoch <- Epoch(pt01EcoG)
fs=1000
sozIndex <- attr(pt01EcoG, "sozIndex")
timeNum <- ncol(epoch)
windowParams<-c(0.25,0.1)

epoch <- Epoch(pt01EcoG)
visuIEEGData(epoch)



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


# mean power analysis over the frequency band over time window (powTimeWindow) and based line time window (baseTimeWindow)
spect<-results[[1]]
timesOnset<-results[[2]]+timeWindow[1]
endBaseIndex<-which.min(abs(timesOnset + 20))

powRange=colSums(spect)
# substract mean baseline t=[-30s:-20s]
mBaseRange<-mean(powRange[1:endBaseIndex])
powRange<-powRange-mBaseRange

