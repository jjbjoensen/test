nm = function(fn, par0 = NULL, lower, upper, N = 1e4, tol = 1e-10) {
  alpha = 1
  gamma = 2
  rho   = 0.5
  sigma = 0.5

  # Create initial simplex
  if (is.null(par0)) {
    n = length(lower)
    simplex = matrix(NA, nrow = n + 1, ncol = n)
    simplex[1, ] = runif(n, lower, upper)
    for (i in 2:(n + 1)) {
      shift = rep(0, n)
      shift[i - 1] = 0.05 * (upper[i - 1] - lower[i - 1])
      simplex[i, ] = simplex[1, ] + shift
    }
  } else {
    n = length(par0)
    simplex = matrix(NA, nrow = n + 1, ncol = n)
    simplex[1, ] = par0
    for (i in 2:(n + 1)) {
      shift = rep(0, n)
      shift[i - 1] = 0.05 * (upper[i - 1] - lower[i - 1])
      simplex[i, ] = par0 + shift
    }
  }

  simplexes = list()
  f_vals = list()

  for (iter in seq_len(N)) {
    # Evaluate
    f = apply(simplex, 1, fn)

    # Sort
    ord = order(f)
    simplex = simplex[ord, , drop = FALSE]
    f = f[ord]

    f_vals[[iter]] = f
    simplexes[[iter]] = simplex

    # Convergence
    if (sd(f) < tol) break

    # Centroid
    centroid = colMeans(simplex[1:n, , drop = FALSE])

    # Reflection
    x_r = centroid + alpha * (centroid - simplex[n + 1, ])
    f_r = fn(x_r)

    if (f_r < f[1]) {
      # Expansion
      x_e = centroid + gamma * (x_r - centroid)
      f_e = fn(x_e)
      if (f_e < f_r) {
        simplex[n + 1, ] = x_e
      } else {
        simplex[n + 1, ] = x_r
      }
    } else if (f_r < f[n]) {
      # Accept reflection
      simplex[n + 1, ] = x_r
    } else {
      # Contraction
      if (f_r < f[n + 1]) {
        x_c = centroid + rho * (x_r - centroid)
      } else {
        x_c = centroid + rho * (simplex[n + 1, ] - centroid)
      }
      f_c = fn(x_c)
      if (f_c < f[n + 1]) {
        simplex[n + 1, ] = x_c
      } else {
        # Shrink towards best
        for (i in 2:(n + 1)) {
          simplex[i, ] = simplex[1, ] + sigma * (simplex[i, ] - simplex[1, ])
        }
      }
    }
  }

  list(simplexes = simplexes, f_vals = f_vals)
}

create_nm_object = function(fn, par0, lower, upper, N, iterations, best_par,
                            best_val, simplex_history, val_history) {
  structure(
    list(
      fn = fn,
      par0 = par0,
      lower = lower,
      upper = upper,
      N = N,
      iterations = iterations,
      best_par = best_par,
      best_val = best_val,
      simplex_history = simplex_history,
      val_history = val_history
    ),
    class = "nm_object"
  )
}


nelder_mead = function(fn, par0 = NULL, lower, upper, N = 1e4, tol = 1e-10) {


  nm_run = nm(fn, par0, lower, upper, N, tol)
  simplex_history = nm_run$simplexes
  val_history = nm_run$f_vals
  iterations = length(val_history)
  if (is.null(par0) == TRUE) {
    initial_pars = "Not provided"
  } else {
    initial_pars = par0
  }

  best_val = nm_run$f_vals[[length(nm_run$f_vals)]][[1]]
  best_par = nm_run$simplexes[[length(nm_run$f_vals)]][1,]


  nm_object = create_nm_object(fn,
                               initial_pars,
                               lower,
                               upper,
                               N,
                               iterations,
                               best_par,
                               best_val,
                               simplex_history,
                               val_history)

  structure(nm_object, class = "nm_object")
}

#' @export
print.nm_object = function(obj) {
  cat("Nelder-Mead Summary:\n")
  cat("Startion parameters:", obj$par0, "\n")
  cat("Lower bound:", obj$lower, "\n")
  cat("Upper bound:", obj$upper, "\n")
  cat("Iterations:", obj$iterations, "\n")
  cat("Best parameters:", obj$best_par, "\n")
  cat("Best function value", obj$best_val, "\n")


  invisible(obj)
}

#' @export
plot.nm_object = function(obj) {
  val_hist = obj$val_history
  val_hist_vec = sapply(val_hist, function(x) x[1])
  print(val_hist_vec)
  point_hist = do.call(rbind, lapply(obj$simplex_history, function(x) x[1, ]))

  d = ncol(point_hist)

  # Plot layout: one for val history + one for each coordinate
  old_par = par(mfrow = c(d + 1, 1), mar = c(4, 4, 2, 1))

  # Plot objective value over iterations
  plot(val_hist_vec, type = "l",
       xlab = "Iteration",
       ylab = "Best Objective Value",
       main = "Objective Value History")

  # Plot each coordinate

  for (j in seq_len(d)) {
    plot(point_hist[, j], type = "l",
         xlab = "Iteration",
         ylab = paste("Coordinate", j),
         main = paste("Evolution of Coordinate", j, "of Best Point"))
  }

  par(old_par)  # Reset layout
}




