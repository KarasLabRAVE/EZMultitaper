standardizeIEEG <- function(data) {
  scaling <- 10^floor(log10(max(data)))
  plotData <- data / scaling
}

#' compute the power analysis over the frequency band using the multitaper method of an epoch list
#'
#' @param epochList list of Matrix or Epoch object. iEEG data matrix or Epoch object. If matrix, the row names are the electrode names and the column names are the time points
#' @param rangeBand Numeric. Frequency band of the multitaper power analysis
#' @param timeBandwidth (numeric): time-half bandwidth product (window duration*half bandwidth of main lobe)
#'                                   (default: 5 Hz*s)
#' @param numTapers (numeric): number of DPSS tapers to use (default: NULL will be computed
#'                                                               as floor(2*time_bandwidth - 1))
#' @param windowParams  (numeric vector): c(window size (seconds), step size (seconds)) (default: 5 1)
#' @param minNfft (numeric) minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
#' @param weighting (char) weighting of tapers ('unity' (default), 'eigen', 'adapt')
#' @param detrendOpt (char): detrend data window ('linear' (default), 'constant', 'off')
#' @param parallel (logical): use parallel processing to speed up calculation (default: FALSE). Note: speedup is faster on
#                             unix-like machines (Mac, Linux) because they allow fork processes while Windows does not.
#         num_workers (numeric): number of cpus/workers to dedicate to parallel processing (default: FALSE). Note: Will
#                                be ignored if parallel is FALSE. If parallel is TRUE and num_workers is false (or if num_workers
#                                exceeds available workers), will default to max number of workers available minus 1.
#' @param numWorkers (logical): number of workers
#'
#' @return A dataframe of power analysis
#' @export
#'
#'@examples
#'# access package Data in OSF
#'dl<-EpochDownloader(progress = FALSE)
#'names(dl)
#'# data download one patient multiple seizure
#'pt01sz1 <- dl$FragilityData_subpt01_1
#'pt01sz1
#'pt01sz2 <- dl$FragilityData_subpt01_2
#'pt01sz2
#'
#'# crop Baseline and signal to analyze
#'useBaseline=TRUE
#'pt01sz1Baseline<-crop(pt01sz1, start=-30, end=-20)
#'pt01sz1Cropped<-crop(pt01sz1, start=0, end=20)
#'pt01sz2Baseline<-crop(pt01sz2, start=-30, end=-20)
#'pt01sz2Cropped<-crop(pt01sz2, start=0, end=20)
#'conditions<-c("30s before sz1","sz1","30s before sz2","sz2")
#'#To visualize patient 01 seizure 1:
#'pt01sozName <- rownames(pt01sz1Cropped)[rowData(pt01sz1Cropped)$soz]
#'pt01Display <- c(pt01sozName, "MLT1", "MLT2", "MLT3", "MLT4")
#'
#'pt01sz1Reordered <- pt01sz1Cropped[pt01Display, ]
#'pt01sz1BSReordered <- pt01sz1Baseline[pt01Display, ]
#'
#'pt01sz2Reordered <- pt01sz2Cropped[pt01Display, ]
#'pt01sz2BSReordered <- pt01sz2Baseline[pt01Display, ]
#'
#'# create Epoch list to follow RAVE/FREEZ framework (have all frequency range, all seizure results ready in UI background)
#'epochList<-list(pt01sz1BSReordered,pt01sz1Reordered,pt01sz2BSReordered,pt01sz2Reordered)
#'
#'names(epochList)<-conditions
generateMultitaper<-function(epochList, rangeBand,
                   timeBandwidth, numTapers, windowParams, minNfft,
                   weighting, detrendOpt, parallel, numWorkers ){

  conditions<-names(epochList)

  df <- data.frame(Conditions = conditions)

  rowIndex <- 1

  fs<-metaData(epochList[[1]])$samplingRate

  for (epoch in epochList){

    timeSeries<-tblData(epoch)

    ne<-nrow(timeSeries)

    spectrogramList <- vector("list",ne)
    cnt <- 1

    for (ie in 1:ne){

      data<-timeSeries[ie,]
      # Compute the multitaper spectrogram
      results <-multitaperSpectrogramR(data=data, fs=fs, frequencyRange=rangeBand,
                                       timeBandwidth=timeBandwidth, numTapers=numTapers, windowParams=windowParams,
                                       minNfft=minNfft, weighting=weighting, detrendOpt=detrendOpt, parallel=parallel, numWorkers=numWorkers)
      spect<-results
      spectrogramList[[cnt]] <- spect
      cnt <- cnt + 1

    }

    #add condition to output
    df$MultitaperData[rowIndex] <- list(spectrogramList)
    rowIndex<-rowIndex+1

  }

  return(df)

}

#' compute the mean power analysis over the frequency band using the multitaper results for a condition
#' for a list of analysis windows
#'
#' @param epochList list of Matrix or Epoch object. iEEG data matrix or Epoch object.
#' @param multitaperResults (numeric)  A dataframe of multitaper power analysis results
#' @param windowParams  (numeric vector): c(window size (seconds), step size (seconds)) (default: 5 1)
#' @param analysisWindows A List of analysis window to generate corresponding heatmap for band mean power analysis
#' @param seizureCondition (char) seizure condition name.
#'
#' @return A list of mean power analysis heatmap
#' @export
#'
#'@examples
#'
#'#Analysis Condition selector
#'seizureCondition="sz1 (2)"
#'useBaseline=TRUE
#'baselineCondition="30s before sz1"
#'baselineDuration<-10

#'#Analysis option
#'#Group1
#'analysisWindows<- vector("list",2)
#'analysisWindows[[1]]$frequencyRange<-c(13,30)
#'analysisWindows[[1]]$timeWindow<-c(0,10)
#'analysisWindows[[2]]$frequencyRange<-c(30,90)
#'analysisWindows[[2]]$timeWindow<-c(0,10)
#'MeanPowBandHeatmaps<-generateHeatmapMeanPowBand(epochList, multitaperResults, analysisWindows, seizureCondition)
generateHeatmapMeanPowBand<-function(epochList, multitaperResults, windowParams, analysisWindows, seizureCondition){

  conditions<-names(epochList)
  condID<-which(conditions==seizureCondition)

  epochSeiz<-epochList[[condID]]
  timeSeries<-tblData(epochSeiz)

  timesSeiz<-as.numeric(colnames(timeSeries))
  elecNames<-rownames(timeSeries)

  # create frequency and time list
  freqList <- list()

  for (i in seq_along(analysisWindows)) {
    freqListName <- paste0("f", i)
    freqList[[freqListName]] <- analysisWindows[[i]]$frequencyRange
  }

  fs<-metaData(epochList[[condID]])$samplingRate
  nel<- nrow(tblData(epochSeiz))
  nt<-ncol(tblData(epochSeiz))
  nwt<-floor((nt/fs-as.numeric(windowParams[1]))/as.numeric(windowParams[2]))+1
  t1<-timesSeiz[1]
  t2<-timesSeiz[length(timesSeiz)]
  heatmapTimes<-round(seq(t1,t2, length.out=nwt),1)

  heatmapFreqList <- vector("list", length(freqList))

  spectrogramList<-multitaperResults$MultitaperData[condID]
  i<-1


for (i in seq_along(freqList)) {

    # extract frequency values
    freqTemp <- freqList[[i]]
    freqStart <- as.numeric(freqTemp[1])
    freqEnd <- as.numeric(freqTemp[2])

    spec<-spectrogramList[[1]][[1]]
    specT <- round(spec[[2]]+t1,1)
    specT <- spec[[2]]+t1
    specF <- spec[[3]]

    heatmap=zeros(nel,nwt)
    freqStartIndex <- which.min(abs(specF - freqStart))
    freqEndIndex  <- which.min(abs(specF - freqEnd))

     for(ie in 1:nel){

        spec<-spectrogramList[[1]][[ie]]
        specR<-spec[[1]]
        spectb <- specR[freqStartIndex:freqEndIndex,]
        heatmapie<-colMeans(spectb)
        heatmap[ie,1:nwt]<-heatmapie[1:nwt]

     }


    meanPowBand<-MeanPowBand(
      pow = heatmap,
      startTimes = specT,
      electrodes = elecNames
    )

    heatmapFreqList[[i]] <- meanPowBand

}

  return(heatmapFreqList)

}


#' compute the mean power analysis over the frequency band using the multitaper method
#'
#' @param epoch Matrix or Epoch object. iEEG data matrix or Epoch object. If matrix, the row names are the electrode names and the column names are the time points
#' @param fs Numeric. frequency of signal iEEG acquisition
#' @param windowParams  Numeric. Window parameters of the multitaper method
#' @param rangeBand Numeric. Frequency band of the multitaper power analysis
#' @param powTimeWindow  Numeric. Time window around seizure onset of the power analysis
#' @param baseTimeWindow Numeric. Time window around seizure onset to subtract the mean baselien signal
#'
#' @return A power analysis matrix
#' @export
#'
#' @examples
#' # access package Data in OSF
#'dl<-EpochDownloader(progress = FALSE)
#'names(dl)
#'# test EZMultitaper with Epoch package
#'# data download one patient one seizure
#'pt01sz1 <- dl$FragilityData_subpt01_1
#'pt01sz1
#'# crop
#'pt01sz1Cropped<-crop(pt01sz1, start=-30, end=20)
#'# select small number of channel
#'pt01sozName <- rownames(pt01sz1Cropped)[rowData(pt01sz1Cropped)$soz]
#'pt01Display <- c(pt01sozName, "MLT1", "MLT2", "MLT3", "MLT4")
#'pt01sz1Reordered <- pt01sz1Cropped[pt01Display, ]
#'betaBand<-c(13,30)
#'rangeBand<-betaBand
#'powTimeWindow<-c(0,20)
#'baseTimeWindow<-c(-30,-20)
#'meta<-metaData(pt01sz1Reordered)
#'fs<-metaData(pt01sz1Reordered)$samplingRate
#'windowParams<-c(0.25,0.1)
#'# compute the mean power analysis over the frequency band (rangeBand) over time window (powTimeWindow) and baselined time window (baseTimeWindow)
#'betaBandPow<-meanPowBaselineBand( epoch=pt01sz1Reordered, fs=fs, windowParams=windowParams, rangeBand=betaBand, powTimeWindow=powTimeWindow, baseTimeWindow=baseTimeWindow)
meanPowBaselineBand <- function(epoch, fs=1000, windowParams, rangeBand, powTimeWindow=c(0,20), baseTimeWindow=c(-30,-20)){


  timeSeries<-tblData(epoch)
  elecNames <- rownames(timeSeries)
  times <- as.numeric(colnames(timeSeries))
  elecNum <- length(elecNames)
  timeNum <- length(times)
  #TODO more option baselining line in mni
  nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1
  #data   <- vector(mode="numeric", length=timeNum)
  data<-timeSeries[1,]
  # Compute the multitaper spectrogram
  results = multitaperSpectrogramR(data=data, fs=fs, windowParams = windowParams, frequencyRange=rangeBand)
  stimes = results[[2]]

  timesOnset<-results[[2]]+times[1]

  startBaseIndex<-which.min(abs(timesOnset - baseTimeWindow[1]))
  endBaseIndex<-which.min(abs(timesOnset - baseTimeWindow[2]))

  elecNum <- nrow(epoch)

  startEpochIndex<-which.min(abs(timesOnset - powTimeWindow[1]))
  endEpochIndex<-which.min(abs(timesOnset - powTimeWindow[2]))

  timesOnset<-timesOnset[startEpochIndex:endEpochIndex]
  powBandMaster=matrix(0,elecNum,length(timesOnset))


  for(ie in 1:elecNum){

    data<-timeSeries[ie,]
    # Compute the multitaper spectrogram
    results <-multitaperSpectrogramR(data=data, fs=fs, windowParams = windowParams, frequencyRange=rangeBand)
    spect<-results[[1]]

    powBandOneElec<-colSums(spect)

    mBasePow<-mean( powBandOneElec[startBaseIndex:endBaseIndex])
    powBandOneElec<- powBandOneElec-mBasePow

    powBandMaster[ie,1:length(timesOnset)]<-powBandOneElec[startEpochIndex:endEpochIndex]

  }

  MeanPowBand(
    pow = powBandMaster,
    startTimes = timesOnset,
    electrodes = elecNames
  )


}


#' Find Seizure Onset Zone
#'
#' The function estimates the seizure onset zone (SOZ). For each row, it calculates the maximum, minimum, or mean of row. The rows with the highest values are considered as the SOZ.
#'
#' @param x MeanPowBand object
#' @param method Character. The method to use to find the onset zone.
#' Must be one of 'max', 'min', or "mean"
#' @param proportion Numeric. The proportion of electrodes to consider as the onset zone.
#' The electrode number will be rounded to the nearest integer.
#' @param ... Additional arguments
#'
#' @return A vector of electrode names, or indices if the electrode names are NULL
#' @export
estimateSOZ <- function(x, method = c("mean", "median", "max", "min"), proportion = 0.1, ...) {
  method <- match.arg(method)

  stopifnot(is(x, "MeanPowBand"))

  powMat <- x$pow
  elCnt <- nrow(powMat)
  nSOZ <- ceiling(elCnt * proportion)
  stopifnot(nSOZ > 0 & nSOZ <= elCnt)

  if (method == "max") {
    stat <- apply(powMat, 1, max)
  } else if (method == "min") {
    stat <- apply(powMat, 1, min)
  } else if (method == "mean") {
    stat <- apply(powMat, 1, mean)
  } else if (method == "median") {
    stat <- apply(powMat, 1, median)
  } else {
    stop("Unknown method: ", method)
  }

  sozIndex <- order(stat, decreasing = TRUE)[seq_len(nSOZ)]
  if (!is.null(x$electrodes)) {
    sozIndex <- x$electrodes[sozIndex]
  }

  sozIndex
}
