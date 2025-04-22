
#' Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
#'
#' @param pow Matrix or mean power object. Either a matrix with row as Electrode names and Column as power index,
#'  or a MeanPowBand object from \code{meanPowBaselineBand}
#'
#' @param sozIndex Integer.  Vector soz electrodes (for good electrodes)
#'
#'
#' @return list of 5 items with quantile matrix, mean and sdv from both electrodes groups
#'
#' @examples
#' data("pt01betaBandPow")
#' data("pt01EcoG")
#' sozIndex <- attr(pt01EcoG, "sozIndex")
#' pt01powstat <- powStat(pow = pt01betaBandPow, sozIndex = sozIndex)
#' @export 
powStat <- function(pow, sozIndex) {
## TODO: support grouped and ungrouped statistics (Not now, but for the future)
    if (is(pow, "MeanPowBand")) pow <- pow$pow
    if (!inherits(pow, "matrix")) stop("pow must be matrix or powility object")
    steps <- ncol(pow)
    sozCID <- which(!(seq_len(nrow(pow)) %in% sozIndex))
    hmapSOZ <- pow[sozIndex, , drop = FALSE]
    hmapREF <- pow[sozCID, , drop = FALSE]
    meanSOZ <- colMeans(hmapSOZ)
    meanRef <- colMeans(hmapREF)
    sdSOZ <- apply(hmapSOZ, 2L, sd)
    sdRef <- apply(hmapREF, 2L, sd)
    Q <- seq(.1, 1, by = .1)
    qmatrix <- rbind(
        apply(hmapSOZ, 2, quantile, Q),
        apply(hmapREF, 2, quantile, Q)
    )
    rowPrefix <- rep(c("SOZ", "REF"), each = 10)
    dimN <- dimnames(qmatrix)
    dimnames(qmatrix) <- list(
        Quantiles = paste0(rowPrefix, dimN[[1L]]),
        Step      = dimN[[2L]]
    )
    PowStat(
        qmatrix   = qmatrix,
        meanSOZ  = meanSOZ,
        meanRef = meanRef,
        sdSOZ    = sdSOZ,
        sdRef   = sdRef
    )
}

