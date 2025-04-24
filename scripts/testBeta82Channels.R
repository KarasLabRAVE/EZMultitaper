
library(R.matlab)
library(readxl)
data <- readMat('data-raw/subpt01sz1_m30sp30s.mat')
pt01EpochRaw <- data$a

## add channel names to the rows
goodChannels <- c(1:4,7:24,26:36,42:43,46:54,56:69,72:95)
sozChannels<-c(33:34,62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01EpochRaw) <- channelNames$name[goodChannels]

sozIndex<-which(goodChannels%in%sozChannels==TRUE)
sozNames<-channelNames$name[sozChannels]

## Add time stamps to the columns
times <- seq(-30, 30, length.out=ncol(pt01EpochRaw))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
colnames(pt01EpochRaw)<-times_with_sign


pt01EcoG<-pt01EpochRaw[,1:50001]
attr(pt01EcoG, "sozIndex") <- sozIndex
attr(pt01EcoG, "sozNames") <- sozNames

fs<-1000
windowParams<-c(0.25,0.1)
betaBand<-c(13,30)
powTimeWindow<-c(0,20)
baseTimeWindow<-c(-30,-20)
epoch <- Epoch(pt01EcoG)
betaBandPow<-meanPowBaselineBand( epoch=epoch, fs=fs, windowParams=windowParams, rangeBand=betaBand, powTimeWindow=powTimeWindow, baseTimeWindow=baseTimeWindow)


## plot the mean power heatmap
plotPowBand<-plotPowHeatmap(pow = betaBandPow, sozIndex = sozIndex)
plotPowBand<-plotPowBand+ggplot2::ggtitle(("Mean beta power heatmap for patient pt01"))
plotPowBand

## plot the mean power quantiles
plotbetaQuantile<-plotPowQuantile(pow = betaBandPow, sozIndex = sozIndex)
plotbetaQuantile<-plotbetaQuantile+ggplot2::ggtitle(("Pooled mean beta power quantiles for patient pt01"))
plotbetaQuantile

## plot the mean power distribution
plotbetaDistr<-plotPowDistribution(pow = betaBandPow, sozIndex = sozIndex)
plotbetaDistr<-plotbetaDistr+ggplot2::ggtitle(("Pooled mean beta power distribution for patient pt01"))
plotbetaDistr