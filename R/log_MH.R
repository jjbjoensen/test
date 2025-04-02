log_mh = function(n, par0, sigmas, ll, lprior, lsampler, lsample_dens) {
  acc = numeric(n)
  m = length(par0)
  pars = matrix(nrow = n, ncol = m)
  par = par0

  for (i in seq_len(n)) {
    park = lsampler(par0, par, sigmas)
    a1 = ll(park) - ll(par)
    a2 = lprior(park, par0, sigmas) - lprior(par, par0, sigmas)
    a3 = lsample_dens(par, par0, park, sigmas) - lsample_dens(park, par0, par, sigmas)

    a = a1 + a2 + a3
    alpha = exp(a)
    u = runif(1)
    if (u <= alpha) {
      acc[[i]] = 1
      pars[i, ] = park
      par = park
    }
    else {
      acc[[i]] = 0
      pars[i, ] = par
    }
  }
  list('pars' = pars, 'acc' = acc)
}

wack = function(n, par0, sigmas, ll, lprior, lsampler, lsample_dens) {
  log_mh_cpp(n, par0, sigmas, ll, lprior, lsampler, lsample_dens)
}
