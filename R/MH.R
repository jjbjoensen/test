mh = function(n, par0, sigmas, likelihood, prior, sampler, sample_dens) {
  acc = numeric(n)
  m = length(par0)
  pars = matrix(nrow = n, ncol = m)
  par = par0

  for (i in seq_len(n)) {
    park = sampler(par0, par, sigmas)
    a1 =  likelihood(park) / likelihood(par)
    a2 = prior(park, par0, sigmas) / prior(par, par0, sigmas)
    a3 = sample_dens(par, par0, park, sigmas) / sample_dens(park, par0, par, sigmas)

    alpha = a1 * a2 * a3
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

