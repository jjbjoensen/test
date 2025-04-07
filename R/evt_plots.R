emp_me = function(data, u) {
  dn = data[data > u]
  1 / length(dn) * sum(data[data > u] - u)
}

me_plot = function(data) {
  data = sort(data)
  data = as.matrix(data)
  N = nrow(data)
  mat = numeric(N)
  for (i in seq_len(nrow(data))) {
    mat[i] = emp_me(data, data[i])
  }
  matrix = cbind(data, mat)
  df = as.data.frame(matrix)
  colnames(df) = c("Value", "ME")
  df
}

std_exp_inverse = function(y) {
  -log(1 - y)
}

qq_against_exponentials = function(data) {
  data_sort = sort(data, decreasing = TRUE)
  n = length(data)
  x = numeric(n)
  y = numeric(n)
  m = n + 1
  for (i in seq_len(n)) {
    x[i] = data_sort[i]
    y[i] = (m - i)
  }
  y = std_exp_inverse(y / m)
  list('x' = x, 'y' = y)
}

# Hill ----

hill_estimator = function(data, k) {
  Xx = sort(data, decreasing = TRUE)
  X = Xx[1:k]
  constant = 1 / k
  inner_sum = numeric(k)
  lXk = log(X[[k]])

  for (i in seq_len(k)) {
    inner_sum[[i]] = log(X[[i]]) - lXk
  }

  (constant * sum(inner_sum))^(-1)

}

hill_plot = function(data) {
  m = length(data)

  y = numeric(m - 1)
  for (i in seq_len(m - 1)) {
    y[[i + 1]] = hill_estimator(data, i + 1)
  }
  x = seq(from = 2, to = m)
  list('x' = x, 'y' = y[2:m])
}

# DEdH -----

H1n = function(data, k) {
  Xx = sort(data, decreasing = TRUE)
  X = Xx[1:(k + 1)]
  constant = 1 / k
  inner_sum = numeric(k)
  lXk = log(X[[k + 1]])

  for (i in seq_len(k)) {
    inner_sum[[i]] = log(X[[i]]) - lXk
  }

  constant * sum(inner_sum)

}


H2n = function(data, k) {
  Xx = sort(data, decreasing = TRUE)
  X = Xx[1:(k + 1)]
  constant = 1 / k
  inner_sum = numeric(k)
  lXk = log(X[[k + 1]])

  for (i in seq_len(k)) {
    inner_sum[[i]] = (log(X[[i]]) - lXk)^2
  }

  constant * sum(inner_sum)
}


DEdH_estimator = function(data, k) {
  H1 = H1n(data, k)
  H2 = H2n(data, k)

  1 + H1 + 1 / 2 * ((H1)^2 / H2 - 1)^(-1)

}

DEdH_plot = function(data) {
  m = length(data)

  y = numeric(m - 1)

  for (i in seq_len(m - 2)) {
    y[[i + 1]] = DEdH_estimator(data, i + 1)
  }
  x = seq(from = 2, to = m)
  list('x' = x, 'y' = y)
}


create_evt_plot_object = function(plot_choice, data, x, y) {
  structure(list(
    plot_choice = plot_choice,
    data = data,
    x = x,
    y = y
  ), class = "evt_plot_object")
}

#' @export
evt_plot = function(plot_choice, data) {
  if (plot_choice == "mean_excess") {
    out = me_plot(data)
    x = out$Value
    y = out$ME
  } else if (plot_choice == "qq_against_exponentials") {
    out = qq_against_exponentials(data)
    x = out$x
    y = out$y
  } else if (plot_choice == "hill") {
    out = hill_plot(data)
    x = out$x
    y = out$y
  } else if (plot_choice == "DEdH") {
    out = DEdH_plot(data)
    x = out$x
    y = out$y
  } else {
    cat("Plot choice not known")
  }



  evt_plot_object = create_evt_plot_object(
    plot_choice = plot_choice,
    data = data,
    x = x,
    y = y
  )

  structure(evt_plot_object, class = "evt_plot_object")
}

#' @export
print.evt_plot_object = function(obj) {
  cat("EVT Plot Summary:\n")
  cat("Chosen plot:", obj$plot_choice, "\n")
  cat("Data summary:", summary(obj$data), "\n")
  cat("X:", head(obj$x), "\n")
  cat("Estimator:", head(obj$y), "\n")


  invisible(obj)
}



#' @export
plot.evt_plot_object = function(obj) {
  if (obj$plot_choice == "mean_excess") {
    plot(
      x = obj$x,
      y = obj$y,
      xlab = "u",
      ylab = "e(u)"
    )

  } else if (obj$plot_choice == "qq_against_exponential") {
    plot(
      x = obj$x,
      y = obj$y,
      xlab = "",
      ylab = ""
    )
  } else if (obj$plot_choice == "Hill") {
    plot(
      x = obj$x,
      y = obj$y,
      xlab = "k",
      ylab = "alpha = 1 / xi",
      type = "l"
    )

  } else if (obj$plot_choice == "DEdH") {
    plot(
      x = obj$x,
      y = obj$y,
      xlab = "k",
      ylab = "xi",
      type = "l"
    )
  } else {
    cat("Plot choice not known")
  }
}
