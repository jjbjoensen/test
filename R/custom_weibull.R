#' @export
rweibull_std = function(n, beta, gamma) {
  samples = numeric(n)
  u = runif(n)
  for (i in seq_len(n)) {
    samples[[i]] = (- log(u[[i]]) / beta)^(1/gamma)
  }
  samples
}

#' @export
dweibull_std = function(x, beta, gamma) {
  beta * gamma * x^(gamma - 1) * exp(- beta * x^gamma)
}

#' @export
pweibull_std = function(q, beta, gamma) {
  1- exp(- beta * q^gamma)
}

#' @export
qweibull_std = function(p, beta, gamma) {
  (- log(1 - p) / beta)^(1 / gamma)
}

#' @export
rweibull_G = function(n, alpha, gamma) {
  samples = numeric(n)
  u = runif(n)
  for (i in seq_len(n)) {
    samples[[i]] = (-log(u[[i]]))^(1/gamma) * alpha
  }
  samples
}

#' @export
dweibull_G = function(x, alpha, gamma) {
  alpha^(-gamma) * gamma * x^(gamma - 1) * exp(- (x / alpha)^gamma)
}

#' @export
pweibull_G = function(q, alpha, gamma) {
  1 - exp(-(q / alpha)^gamma)
}

#' @export
qweibull_G = function(p, beta, gamma) {
  (- log(1 - p))^(1/gamma) / alpha
}


