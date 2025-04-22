library(matrixStats)
data("pt01EcoG")
timeWindow <- c(-30, 20)
epoch <- Epoch(pt01EcoG)
fs=1000
sozIndex <- attr(pt01EcoG, "sozIndex")
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
powTimeWindow<-c(0,20)
baseTimeWindow<-c(-30,-20)

# compute the mean power analysis over the frequency band (rangeBand) over time window (powTimeWindow) and baselined time window (baseTimeWindow)
betaBandPow<-meanPowBaselineBand( epoch=epoch, fs=fs, windowParams=windowParams, rangeBand=betaBand, powTimeWindow=powTimeWindow, baseTimeWindow=baseTimeWindow)

plotPowBand<-plotPowHeatmap(pow=betaBandPow,sozIndex=sozIndex)
plotPowBand


############### test example 1

data("pt01EcoG")
fs<-1000
windowParams<-c(0.25,0.1)
betaBand<-c(13,30)
powTimeWindow<-c(0,20)
baseTimeWindow<-c(-30,-20)
epoch <- Epoch(pt01EcoG)
betaBandPow<-meanPowBaselineBand( epoch=epoch, fs=fs, windowParams=windowParams, rangeBand=betaBand, powTimeWindow=powTimeWindow, baseTimeWindow=baseTimeWindow)
 
############### test example 2
data("pt01EcoG")

## sozIndex is the index of the electrodes we assume are in the SOZ
sozIndex <- attr(pt01EcoG, "sozIndex")
## precomputed MeanPowBand object
data("pt01betaBandPow")

## plot the mean power heatmap
plotPowBand<-plotPowHeatmap(pow = pt01betaBandPow, sozIndex = sozIndex)
plotPowBand<-plotPowBand+ggplot2::ggtitle(("Mean beta power heatmap for patient pt01"))
plotPowBand

## plot the mean power quantiles
plotbetaQuantile<-plotPowQuantile(pow = pt01betaBandPow, sozIndex = sozIndex)
plotbetaQuantile<-plotbetaQuantile+ggplot2::ggtitle(("Pooled mean beta power quantiles for patient pt01"))
plotbetaQuantile

## plot the mean power distribution
plotbetaDistr<-plotPowDistribution(pow = pt01betaBandPow, sozIndex = sozIndex)
plotbetaDistr<-plotbetaDistr+ggplot2::ggtitle(("Pooled mean beta power distribution for patient pt01"))
plotbetaDistr

############# test example 3
data("pt01betaBandPow")
data("pt01EcoG")
sozIndex <- attr(pt01EcoG, "sozIndex")
pt01PowStat <- powStat(pow = pt01betaBandPow, sozIndex = sozIndex)
