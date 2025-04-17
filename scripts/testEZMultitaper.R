data("pt01EcoG")
timeWindow <- c(-30, 20)
epoch <- Epoch(pt01EcoG)
fs=1000
sozIndex <- attr(pt01EcoG, "sozIndex")
timeNum <- ncol(epoch)
windowParams = c(1, 0.2) 

nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1
data   <- vector(mode="numeric", length=timeNum)
data[1:timeNum]<-epoch$data[sozIndex[1],1:timeNum]
# Compute the multitaper spectrogram
results = multitaperSpectrogramR(data=data, fs=fs, windowParams=c(0.25,0.1))
## Visualize a subject of electrodes
#display <- c(sozIndex, 77:80)
stimes = results[[2]]
epoch <- Epoch(pt01EcoG)
visuIEEGData(epoch)

timesOnset<-stimes+timeWindow[1]

nWindow<-length(stimes)

baselinetime<-timesOnset[1:98]
spectrogram<-nanpow2db(results[[1]])

baselinespec<-spectrogram[,1:98]
meanbaseline<-rowMeans(baselinespec)
spectmmean<-spectrogram[1:nrow(spectrogram),]-meanbaseline[1:nrow(spectrogram)]

timesOnset<-timesOnset[99:ncol(spectrogram)]
spectm20sp20s<-spectmmean[,99:ncol(spectrogram)]

spectrogram<-spectm20sp20s

colnames(spectrogram)<-timesOnset
rownames(spectrogram)<-ceiling(results[[3]])

flipspect<-spectrogram[nrow(spectrogram):1,]

dfspec<-data.frame(flipspect)

xlabel<-"Time (s)"
ylabel<-"Frequency (Hz)"

plotf<-makeHeatMap(df=dfspec)+
  ggplot2::labs(x = xlabel, y=ylabel, size = 2) 
plotf


