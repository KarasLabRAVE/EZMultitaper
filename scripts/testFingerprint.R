library(matrixStats)
data("pt01EcoG")
timeWindow <- c(-30, 20)
epoch <- Epoch(pt01EcoG)
fs=1000
sozIndex <- attr(pt01EcoG, "sozIndex")
timeNum <- ncol(epoch)
windowParams<-c(0.25,0.1)
frequencyRange<-c(0.5,200)

nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1
data   <- vector(mode="numeric", length=timeNum)
data[1:timeNum]<-epoch$data[sozIndex[1],1:timeNum]
# Compute the multitaper spectrogram
results = multitaperSpectrogramR(data=data, fs=fs, windowParams = windowParams, frequencyRange=frequencyRange)
## Visualize a subject of electrodes
#display <- c(sozIndex, 77:80)
stimes = results[[2]]
epoch <- Epoch(pt01EcoG)
visuIEEGData(epoch)

timesOnset<-stimes+timeWindow[1]

nWindow<-length(stimes)

nBase<-ceiling(nWindow/5)

baselinetime<-timesOnset[1:nBase]
spectrogram<-nanpow2db(results[[1]])
#spectrogram<-results[[1]]

baselinespec<-spectrogram[,1:nBase]
meanbaseline<-rowMeans(baselinespec)
stdbaseline<-rowSds(baselinespec)
spectbaseline<-matrix(0,nrow(spectrogram),ncol(spectrogram))
for(i in 1:nrow(spectrogram)){
  #spectbaseline[i,]<-(spectrogram[i,]-meanbaseline[i])/stdbaseline[i]
  spectbaseline[i,]<-(spectrogram[i,]-meanbaseline[i])
  
}

#spectbaseline<-matrix_threshold(spectbaseline,minval=0)

nBase<=nBase+1
timesOnset<-timesOnset[nBase:ncol(spectrogram)]
#spectm20sp20s<-spectrogram[,nBase:ncol(spectrogram)]
spectm20sp20s<-spectbaseline[,nBase:ncol(spectrogram)]
spectrogram<-spectm20sp20s

colnames(spectrogram)<-ceiling(timesOnset)
rownames(spectrogram)<-ceiling(results[[3]])

flipspect<-spectrogram[nrow(spectrogram):1,]

dfspec<-data.frame(flipspect)

xlabel<-"Time (s)"
ylabel<-"Frequency (Hz)"

plotf<-makeHeatMapDiscretexy(df=dfspec)+
  ggplot2::labs(x = xlabel, y=ylabel, size = 2) 
plotf