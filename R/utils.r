
# ---- Helpers ----
`%||%` <- function(a, b) if (!is.null(a)) a else b
label01 <- function(x) ifelse(toupper(trimws(x)) == "F", 1L, 0L)


to_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  dimnames(x) <- NULL
  x
}


#' Check and keep valid index only
#'
#' @param indices Numeric or character index to check
#' @param names Character. All names corresponding to the indices
checkIndex <- function(indices, names) {
  if (length(names) == 0) {
    return()
  }
  if (length(indices) == 0) {
    return()
  }
  if (is(indices, "numeric")) {
    allIndices <- seq_along(names)
    diffIndices <- setdiff(indices, allIndices)
    indicesFiltered <- indices[!indices %in% diffIndices]
    result <- indicesFiltered
  } else {
    diffIndices <- setdiff(indices, names)
    indicesFiltered <- indices[!indices %in% diffIndices]
    result <- which(names %in% indicesFiltered)
  }
  if (length(diffIndices)) {
    indicesMissing <- paste(diffIndices, collapse = ", ")
    indicesExist <- paste(indicesFiltered, collapse = ", ")
    warning(
      glue("Indices {indicesMissing} are out of range. I will keep the valid values {indicesExist}.")
    )
  }
  result
}


#pool score
.clamp01 <- function(p) pmin(1 - 1e-6, pmax(1e-6, p))
.logit   <- function(p) log(.clamp01(p) / (1 - .clamp01(p)))
.invlogit<- function(x) 1 / (1 + exp(-x))


