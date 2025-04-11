create_mh_object = function(par0, mh_choice, params_history, acc_history, moving_averages) {
  structure(
    list(
      par0 = par0,
      sampler_choice = mh_choice,
      params_history = params_history,
      acc_history = acc_history,
      moving_averages = moving_averages
    ),
    class = "mh_object"
  )
}

mh_validator = function(n, par0, sigmas, likelihood, lprior, sampler, lsample_dens) {
  stopifnot("Provided n cannot be converted to integer"=is.integer(as.integer(n)))
  stopifnot("Provided par0 is not a vector"=is.vector(par0))
  stopifnot("Provided (log-)likelihood is not a function"=is.function(likelihood))
  stopifnot("Provided (log-)prior is not a function"=is.function(lprior))
  stopifnot("Provided sampler is not a function"=is.function(sampler))
  stopifnot("Provided sample (log-)density is not a function"=is.function(lsample_dens))
}

#' @export
metropolis_hastings = function(n, par0, sigmas, likelihood, prior, sampler,
                               sample_dens, log = FALSE) {

  mh_validator(n, par0, sigmas, likelihood, prior, sampler, sample_dens)
  if (log == TRUE) {
    mh_choice = log_mh_cpp
    sampler_choice = "log"
  } else {
    mh_choice = mh_cpp
    sampler_choice = "standard"
  }

  mh_run = mh_choice(n, par0, sigmas, likelihood, prior, sampler, sample_dens)
  params_history = mh_run$parameters
  acc_history = mh_run$acceptance
  m = length(par0)
  avg = matrix(nrow = n, ncol = m)
  for (i in seq_len(m)) {
    avg[, i] = moving_average(params_history[, i])
  }
  mh_object = create_mh_object(
    par0 = par0,
    mh_choice = sampler_choice,
    params_history = params_history,
    acc_history = acc_history,
    moving_averages = avg
  )

  structure(mh_object, class = c("mcmc_object", "mh_object"))
}

#' @export
print.mh_object = function(obj) {
  cat("MH Summary:\n")
  cat("MH Type:", obj$sampler_choice[1], "\n")
  cat("Initial Parameters:", obj$par0, "\n")
  cat("Acceptance ratio:", mean(obj$acc_history), "\n")
  cat("Iterations:", length(obj$acc_history), "\n")


  invisible(obj)
}

moving_average = function(x) {
  n = length(x)
  ones = rep(1, n)
  ma = cumsum(x) / cumsum(ones)
  ma
}

#' @export
plot.mh_object = function(obj) {
  n = length(obj$par0)
  params_history = obj$params_history
  plots = list()

  for (i in seq_len(n)) {
    avg = moving_average(params_history[, i])

    plots[[i]] = plot(avg, type = "l",
                      main = paste("Parameter", i, "Moving Average"),
                      xlab = "Iteration", ylab = paste("Parameter", i))
  }

}


#' @export
summarise.mh_object = function(obj, burn_in) {
  n = length(obj$par0)
  map = numeric(n)
  bayes = numeric(n)
  acc_ratio = numeric(n)
  params_history = obj$params_history
  acc_history = obj$acc_history
  acc_ratio = mean(acc_history)
  m = length(params_history)

  for (i in seq_len(n)) {
    map[[i]] = estimate_mode(params_history[, i][burn_in:m])
    bayes[[i]] = mean(params_history[, i][burn_in:m])
  }

  list('MAP' = map, 'Bayes' = bayes, 'Acc_ratio' = acc_ratio)

}

#' @export
optimize_mh = function(n_range = 1e4, n_optim = 1e3, optim_maxit = 50,
                       sigmas0, par0, likelihood, prior, sampler,
                       prior_optim = F, sample_dens, log = FALSE) {

  range_tries = c(1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
  m = length(range_tries)
  acc_ratios = numeric(m)
  k = length(par0)

  for (i in seq_len(m)) {
    if (prior_optim == FALSE) {
      sigmas = c(rep(range_tries[[i]],k),tail(sigmas0, k))
    } else {
      sigmas = rep(range_tries,length(sigmas0))
    }

    cat("sigmas", sigmas, "\n")
    tests = metropolis_hastings(n_range, par0, sigmas, likelihood, prior,
                                sampler, sample_dens, log)
    acc_ratios[[i]] = sum(tests$acc_history)
  }

  if (prior_optim == FALSE) {
    sigma_range = c(rep(range_tries[which.max(acc_ratios)],k),tail(sigmas0, k))
  } else {
    sigma_range = rep(range_tries[which.max(acc_ratios)],length(sigmas0))
  }

  fn = function(sigmas) {
    tests =  metropolis_hastings(n_optim, par0, sigmas, likelihood, prior,
                                 sampler, sample_dens, log)
    m = sum(tests$acc_history)

    -1 * m
  }

  best_sigma = optim(sigma_range, fn, method = "Nelder-Mead",
                     control = list(optim_maxit))

  list('best_sigmas'  = best_sigma$par, 'acc_ratio' = best_sigma$value / n_optim)

}


