#' @export
rpareto = function(n, alpha, theta) {
  samples = numeric(n)
  u = runif(n)
  for (i in seq_len(n)) {
    samples[[i]] = theta * ((1 - u[[i]])^(-1/alpha) - 1)
  }
  samples
}

#' @export
dpareto = function(x, alpha, theta) {
  alpha * theta^(alpha) / ((theta + x)^(1+alpha))
}

#' @export
ppareto = function(q, alpha, theta) {
  1 - (theta / (theta + q))^alpha
}

#' @export
qpareto = function(p, alpha, theta) {
  theta * ((1 - p)^(-1/alpha) - 1)
}

