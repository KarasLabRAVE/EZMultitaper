library(EZMultitaper)
library(Epoch)
library(ggplot2)
library(ggtext)
library(gsignal)

# access package Data in OSF
dl<-EpochDownloader(progress = FALSE)
names(dl)

# test EZMultitaper with Epoch package
# data download one patient one seizure
# adapt for RAVE use with voltage list
pt01sz1 <- dl$FragilityData_subpt01_1
pt01sz2 <- dl$FragilityData_subpt01_2

#pt01sz1ieeg<-tblData(pt01sz1)

# crop Baseline and signal to analyze
pt01sz1Baseline<-crop(pt01sz1, start=-30, end=-20)
pt01sz1Cropped<-crop(pt01sz1, start=-10, end=20)
pt01sz2Baseline<-crop(pt01sz2, start=-30, end=-20)
pt01sz2Cropped<-crop(pt01sz2, start=-10, end=20)

gainInput=0.5
conditions<-c("30s before sz1","sz1 (2)","30s before sz2 (3)","sz2 (4)")

# 120925 Load and build Epoch Package
# from https://github.com/KarasLabRAVE/Epoch/tree/devel-cecile
# Change gain
visualSOZ <- function(epoch, sozNames, gainInput) {


  p <- plot(epoch,gain=gainInput)

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
pt01sz1BSReordered <- pt01sz1Baseline[pt01Display, ]

pt01sz2Reordered <- pt01sz2Cropped[pt01Display, ]
pt01sz2BSReordered <- pt01sz2Baseline[pt01Display, ]

# create Epoch list to follow RAVE/FREEZ framework (have all frequency range, all seizure results ready in UI background)
epochList<-list(pt01sz1BSReordered,pt01sz1Reordered,pt01sz2BSReordered,pt01sz2Reordered)

names(epochList)<-conditions

names(epochList)
epoch<-epochList[[4]]
visualSOZ(epoch=epoch,pt01sozName,gainInput)

#Set default spectrogram parameters
multitaperBand<-c(0.5,250)
minNfft<-0.0
numWorkers<-31
parallel<-FALSE
detrendOpt<-'linear'
weighting<-'unity'
windowParams<-c(2.5,0.5)
timeBandwidth<-3.0
numTapers<-5
meta<-metaData(pt01sz1Reordered)
fs<-metaData(pt01sz1Reordered)$samplingRate
windowParams<-c(0.25,0.1)
timeSeries<-tblData(epoch)

# compute the multitaper results for each electrodes of each epoch
multitaperResults<-generateMultitaper(epochList, rangeBand=multitaperBand,
                   timeBandwidth, numTapers, windowParams, minNfft,
                   weighting, detrendOpt, parallel, numWorkers )

#Analysis Condition selector
seizureCondition="sz1 (2)"
useBaseline=TRUE
baselineCondition="30s before sz1"
baselineDuration<-10

#Analysis option
#Group1
analysisWindows<- vector("list",2)
analysisWindows[[1]]$frequencyRange<-c(13,30)
analysisWindows[[1]]$timeWindow<-c(0,10)
analysisWindows[[2]]$frequencyRange<-c(30,90)
analysisWindows[[2]]$timeWindow<-c(0,10)

MeanPowBandHeatmaps<-generateHeatmapMeanPowBand(epochList, multitaperResults, windowParams, analysisWindows, seizureCondition)


# power heatmap with the same display options as the previous voltage plot.
# Looking at both plots allows to check correlation between soz patterns
powHeatmap <- function(heatmap, sozNames) {
  startTimes <- heatmap$startTimes

  indexsz <- which(abs(startTimes)<=0.01)
  elecColor <- rep("black", length(heatmap$electrodes))
  elecColor[heatmap$electrodes%in% sozNames] <- 'red'
  elecColor <- rev(elecColor)

  plot(heatmap)+
    # geom_vline(xintercept = indexsz, color = "black", linetype = "dashed", linewidth = 1) +
    theme(
      axis.text.y = element_markdown(colour = elecColor)
    )
}

powHeatmap(heatmap=MeanPowBandHeatmaps[[1]], sozNames=pt01sozName)


## plot the mean band power distribution
plotPowDistribution(pow = MeanPowBandHeatmaps[[1]], groupIndex = pt01sozName, bandType="SEM", rollingWindow = 1)

plotPowQuantile(pow = MeanPowBandHeatmaps[[1]], groupIndex = pt01sozName)
