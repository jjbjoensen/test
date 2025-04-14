#' @export
estimate_mode = function(x) {
  d = density(na.omit(x))
  d$x[which.max(d$y)]
}
