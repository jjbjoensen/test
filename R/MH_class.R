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
  stopifnot("Provided (log-)likelihood is not a function"=is.function(ll))
  stopifnot("Provided (log-)prior is not a function"=is.function(lprior))
  stopifnot("Provided sampler is not a function"=is.function(sampler))
  stopifnot("Provided sample (log-)density is not a function"=is.function(lsample_dens))
}

#' @export
metropolis_hastings = function(n, par0, sigmas, likelihood, prior, sampler,
                               sample_dens, log = FALSE) {

  mh_validator(n, par0, sigmas, ll, prior, sampler, sample_dens)
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


