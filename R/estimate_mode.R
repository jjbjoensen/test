#' @export
estimate_mode = function(x) {
  d = density(x)
  d$x[which.max(d$y)]
}
