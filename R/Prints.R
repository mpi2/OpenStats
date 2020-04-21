print.NULL <- function(x, ...) {
  message0("No print available for NULL objects")
}
print.OpenStatsFE <- function(x, format = "rst", ...) {
  r <- summary.OpenStatsFE(object = x, format = format, ...)
  invisible(r)
}
print.OpenStatsRR <- function(x, format = "rst", ...) {
  r <- summary.OpenStatsRR(object = x, format = format, ...)
  invisible(r)
}
print.OpenStatsMM <- function(x, format = "rst", ...) {
  r <- summary.OpenStatsMM(object = x, format = format, ...)
  invisible(r)
}
