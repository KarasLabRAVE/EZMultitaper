
# ---- Helpers ----
`%||%` <- function(a, b) if (!is.null(a)) a else b
label01 <- function(x) ifelse(toupper(trimws(x)) == "F", 1L, 0L)


to_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  dimnames(x) <- NULL
  x
}

#pool score
.clamp01 <- function(p) pmin(1 - 1e-6, pmax(1e-6, p))
.logit   <- function(p) log(.clamp01(p) / (1 - .clamp01(p)))
.invlogit<- function(x) 1 / (1 + exp(-x))
