# test EZMultitaper with Epoch package
# data download one patient one seizure
dl <- EpochDownloader()
pt01sz1<-dl$Retrostudy_subpt01_1

#pt01sz1ieeg<-tblData(pt01sz1)

# crop
pt01sz1m30p20s<-crop(pt01sz1, start=-30, end=20)
sozIndex<-which(pt01sz1m30p20s@rowData$soz==TRUE)

display <- c(sozIndex, 77:80)
# subset on 14 electrodes
pt01sz1m30p20s14e<-pt01sz1m30p20s[display,]
sozIndex14e<-which(pt01sz1m30p20s14e@rowData$soz==TRUE)

# in this package visuiEEG uses the epoch package
visuIEEGData(epoch=pt01sz1m30p20s14e)

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

fs=1000
windowParams<-c(0.25,0.1)


# compute the mean power analysis over the frequency band (rangeBand) over time window (powTimeWindow) and baselined time window (baseTimeWindow)
betaBandPow<-meanPowBaselineBand( epoch=pt01sz1m30p20s14e, fs=fs, windowParams=windowParams, rangeBand=betaBand, powTimeWindow=powTimeWindow, baseTimeWindow=baseTimeWindow)




## plot the mean power heatmap
plotPowBand<-plotPowHeatmap(pow = betaBandPow, sozIndex = sozIndex14e)
plotPowBand<-plotPowBand+ggplot2::ggtitle(("Mean beta power heatmap for patient pt01"))
plotPowBand

## plot the mean power quantiles
plotbetaQuantile<-plotPowQuantile(pow = betaBandPow, sozIndex = sozIndex14e)
plotbetaQuantile<-plotbetaQuantile+ggplot2::ggtitle(("Pooled mean beta power quantiles for patient pt01"))
plotbetaQuantile

## plot the mean power distribution
plotbetaDistr<-plotPowDistribution(pow = betaBandPow, sozIndex = sozIndex14e)
plotbetaDistr<-plotbetaDistr+ggplot2::ggtitle(("Pooled mean beta power distribution for patient pt01"))
plotbetaDistr

pt01PowStat <- powStat(pow = betaBandPow, sozIndex = sozIndex14e)
