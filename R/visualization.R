# A plot function that takes a data frame and returns a heatmap plot
makeHeatMap <- function(df, xTicksNum = 10, maxLabels = Inf){
  xLabels <- colnames(df)
  yLabels <- rownames(df)

  if(is.null(xLabels)){
    xLabels <- seq_len(ncol(df))
  }
  if(is.null(yLabels)){
    yLabels <- seq_len(nrow(df))
  }

  df$y <- yLabels
  df_long <- reshape2::melt(df, id.vars = "y", variable.name = "x", value.name = "value")
  colnames(df_long) <- c("y", "x", "value")

  ## sort df_long by rownames
  df_long$x <- factor(df_long$x, levels = xLabels)
  df_long$y <- factor(df_long$y, levels = rev(yLabels))
  ## show 10 time points on x-axis at most
  if (length(xLabels) > xTicksNum){
    step <- ceiling(length(xLabels) / xTicksNum)
    breaksIdx <- seq(1, length(xLabels), by = step)
    #breaks <- as.character(round(as.numeric(xLabels[breaksIdx]),digits=1))
    breaks <- xLabels[breaksIdx]
  } else {
    breaks <- xLabels
  }

  ## limit the number of labels on y-axis
  yLabelsForDisplay <- rev(yLabels)  # Match the reversed factor levels
  if (length(yLabelsForDisplay) > maxLabels) {
    by_num <- ceiling(length(yLabelsForDisplay)/maxLabels)
    label_idx <- seq(length(yLabelsForDisplay), 1, by=-by_num)
    yLabelsForDisplay[-label_idx] <- ""
  }

  ggplot(df_long) +
    geom_tile(aes(x = .data$x, y = .data$y, fill = .data$value)) +
    scale_x_discrete(labels = breaks, breaks = breaks) +
    scale_y_discrete(labels = yLabelsForDisplay, breaks = rev(yLabels)) +
    theme(plot.title = element_markdown(hjust = 0.5)) +
    viridis::scale_fill_viridis(option = "turbo") +
    theme_minimal()
}



#' Visualization functions (raw signal, mean power band matrix)
#'
#' @description `plot`: plot mean power band with electrodes marked as soz colored
#'
#' @param x MeanPowBand object from \code{meanPowBaselineBand}
#' @param y Not used (for S4 method compatibility)
#' @param pow meanPowBaselineBand object from \code{meanPowBaselineBand}
#' @param x.lab.size Numeric. Size of x-axis labels. Default is 4.
#' @param y.lab.size Numeric. Size of y-axis labels. Default is 10
#' @inheritParams powStat
#' @param maxLabels Integer. Maximum number of labels to show on y-axis. Default is 50. The actual number of labels may be less than this value if there are too many electrodes.
#'
#' @return A ggplot object
#'
#' @examples
#'
#' data("pt01EcoG")
#' @rdname plotPow
#' @export
setMethod("plot", signature(x = "MeanPowBand", y = "missing"),
          function(x, y,
                   groupIndex = NULL,
                   maxLabels = 50,
                   ranked = FALSE,
                   x.lab.size = 10,
                   y.lab.size = 10) {
            powMat <- x$pow

            elecNum <- nrow(powMat)
            windowNum <- ncol(powMat)

            elecNames <- x$electrodes
            groupIndex <- checkIndex(groupIndex, elecNames)

            group1 <- groupIndex
            group2 <- setdiff(seq_len(elecNum), groupIndex)

            elecColor <- rep("blue", elecNum)
            elecColor[seq_along(group2)] <- "black"

            startTime <- x$startTimes
            if (is.null(startTime)) {
              xlabel <- "Time Index"
              stimes <- seq_len(windowNum)
            } else {
              xlabel <- "Time (s)"
              stimes <- round(startTime,digits=1)
            }

            rownames(powMat) <- x$electrodes
            colnames(powMat) <- stimes

            ## prepare the data.frame for visualization
            allIndex <- c(group1, group2)
            df <- as.data.frame(powMat[allIndex, ])

            makeHeatMap(df, maxLabels = maxLabels) +
              labs(x = xlabel, y = "Electrode") +
              theme(
                axis.text.x = element_text(size = x.lab.size),
                axis.text.y = element_markdown(size = y.lab.size, colour = elecColor), # Adjust depending on electrodes
              )
          }
)






#' @description `plotPowQuantile`: Plot mean power time quantiles for two electrodes group marked as SOZ and reference
#'
#' @rdname plotPowHeatmap
#' @examples
#' ## plot the mean power quantiles
#' plotbetaQuantile<-plotPowQuantile(pow = pt01betaBandPow, sozIndex = sozIndex)
#' plotbetaQuantile
#' @export
plotPowQuantile <- function(pow, sozIndex = NULL) {
  sozIndex <- checkIndex(sozIndex, pow$electrodes)
  if (is.null(sozIndex)) {
    sozIndex <- estimateSOZ(pow)
  }
  windowNum <- ncol(pow)

  stat <- powStat(pow, sozIndex)
  qmatrix <- as.data.frame(stat$qmatrix)

  startTimes <- pow$startTimes
  if (is.null(startTimes)) {
    xlabel <- "Time Index"
    timeTicks <- seq_len(windowNum)
  } else {
    xlabel <- "Time (s)"
    timeTicks <- startTimes
  }

  colnames(qmatrix) <- timeTicks

  makeHeatMap(qmatrix)+
    ggplot2::labs(x = xlabel, y = "Quantiles", size = 2) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 4), # Adjust depending on electrodes
    )
}


#' @description `plotPowDistribution`: Plot mean power time distribution for two electrodes group marked as SOZ and reference
#'
#' @rdname plotPowHeatmap
#' @examples
#' ## plot the mean power distribution
#' plotBetaDistr<-plotPowDistribution(pow = pt01betaBandPow, sozIndex = sozIndex)
#' plotBetaDistr
#' @export
plotPowDistribution <- function(pow, sozIndex = NULL) {
  if (is.null(sozIndex)) {
    sozIndex <- estimateSOZ(pow)
  }

  sozIndex <- checkIndex(sozIndex, pow$electrodes)

  powMat <- pow$pow
  windowNum <- ncol(powMat)

  SOZMat <- powMat[sozIndex, , drop = FALSE]
  RefMat <- powMat[-sozIndex, , drop = FALSE]

  meanSOZ <- apply(SOZMat, 2, mean, na.rm = TRUE)
  semSOZ <- apply(SOZMat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

  meanRef <- apply(RefMat, 2, mean, na.rm = TRUE)
  semRef <- apply(RefMat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

  startTimes <- pow$startTimes
  if (is.null(startTimes)) {
    xlabel <- "Time Index"
    timeTicks <- seq_len(windowNum)
  } else {
    xlabel <- "Time (s)"
    timeTicks <- startTimes
  }

  upperSOZ <- meanSOZ + semSOZ
  lowerSOZ <- meanSOZ - semSOZ
  upperRef <- meanRef + semRef
  lowerRef <- meanRef - semRef

  plotData <- data.frame(
    timeTicks = timeTicks,
    meanSOZ = meanSOZ,
    upperSOZ = upperSOZ,
    lowerSOZ = lowerSOZ,
    meanRef = meanRef,
    upperRef = upperRef,
    lowerRef = lowerRef
  )

  colors <- c("SOZ +/- sem" = "red", "SOZc +/- sem" = "black")
  ggplot2::ggplot(plotData, ggplot2::aes(x = .data$timeTicks)) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab("Mean Power") +
    ggplot2::geom_line(ggplot2::aes(y = .data$meanSOZ, color = "SOZ +/- sem")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$upperSOZ), color = "red", linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(y = .data$lowerSOZ), color = "red", linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(y = .data$meanRef, color = "SOZc +/- sem")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$upperRef), color = "black", linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(y = .data$lowerRef), color = "black", linetype = "dotted") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowerSOZ, ymax = .data$upperSOZ), fill = "red", alpha = 0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowerRef, ymax = .data$upperRef), fill = "black", alpha = 0.5) +
    ggplot2::scale_color_manual(name = "Electrode groups", values = c(colors))
}
