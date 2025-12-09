library(Epoch)
library(ggplot2)
library(ggtext)
library(gsignal)

# access package Data in OSF
dl<-EpochDownloader(progress = FALSE)
names(dl)

# test EZMultitaper with Epoch package
# data download one patient one seizure
pt01sz1 <- dl$FragilityData_subpt01_1
pt01sz1

#pt01sz1ieeg<-tblData(pt01sz1)

# crop
pt01sz1Cropped<-crop(pt01sz1, start=-30, end=20)

# Change gain
visualSOZ <- function(epoch, sozNames) {
  p <- plot(epoch)

  elecColor <- rep("black", nrow(epoch))
  elecColor[rownames(epoch) %in% sozNames] <- 'red'
  elecColor <- rev(elecColor) # match the electrode order in the plot

  p +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1)+
    theme(axis.text.y = element_markdown(colour = elecColor))
}

#To visualize patient 01 seizure 1:
pt01sozName <- rownames(pt01sz1Cropped)[rowData(pt01sz1Cropped)$soz]
pt01Display <- c(pt01sozName, "MLT1", "MLT2", "MLT3", "MLT4")
pt01sz1Reordered <- pt01sz1Cropped[pt01Display, ]

visualSOZ(pt01sz1Reordered, pt01sozName)

#Set spectrogram parameters
#define frequency bands
deltaBand<-c(0.5,4)
thetaBand<-c(4,8)
alphaBand<-c(8,13)
betaBand<-c(13,30)
gammaBand<-c(30,90)
highGammaBand<-c(90,150)

rangeBand<-betaBand
powTimeWindow<-c(0,20)
baseTimeWindow<-c(-30,-20)

meta<-metaData(pt01sz1Reordered)
fs<-metaData(pt01sz1Reordered)$samplingRate
windowParams<-c(0.25,0.1)


# compute the mean power analysis over the frequency band (rangeBand) over time window (powTimeWindow) and baselined time window (baseTimeWindow)
betaBandPow<-meanPowBaselineBand( epoch=pt01sz1Reordered, fs=fs, windowParams=windowParams, rangeBand=betaBand, powTimeWindow=powTimeWindow, baseTimeWindow=baseTimeWindow)

## plot the mean power heatmap
plotPowBand<-plotPowHeatmap(pow = betaBandPow)
plotPowBand<-plotPowBand+ggplot2::ggtitle(("Mean beta power heatmap for patient pt01"))
plotPowBand

# ## plot the mean power quantiles
# plotbetaQuantile<-plotPowQuantile(pow = betaBandPow)
# plotbetaQuantile<-plotbetaQuantile+ggplot2::ggtitle(("Pooled mean beta power quantiles for patient pt01"))
# plotbetaQuantile
#
# ## plot the mean power distribution
# plotbetaDistr<-plotPowDistribution(pow = betaBandPow, sozIndex = sozIndex14e)
# plotbetaDistr<-plotbetaDistr+ggplot2::ggtitle(("Pooled mean beta power distribution for patient pt01"))
# plotbetaDistr
#
# pt01PowStat <- powStat(pow = betaBandPow, sozIndex = sozIndex14e)
