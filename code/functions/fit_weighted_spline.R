fit_weighted_spline <- function(x, params, weights, smoothing_vals) {
  require(furrr, quietly = TRUE)
  require(future, quietly = TRUE)
  require(magrittr, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(purrr, quietly = TRUE)

  coefs <- list()
  d <- ncol(params)
  smoothing_error <- vector(mode = "numeric", length = length(smoothing_vals))
  lambda <- smoothing_vals

  train_param_mat <- cbind(rep(1, nrow(params)), params)
  train_weights <- weights
  train_E <- calcE(params, 4 - d)

  best_error <- 1e10
  best_spline <- NULL
  opt_run <- 1

  for (smoothing_idx in seq_along(smoothing_vals)) {
    spline_est <- solve_weighted_spline_hat(
      train_E,
      train_weights,
      train_param_mat,
      x,
      lambda[smoothing_idx],
      ncol(params),
      ncol(x)
    )

    fitted_values <- spline_est$hat %*% x

    RSS <- apply(x - fitted_values, 1, norm_euclidean)^2 |>
      sum()

    gcv_denom <- nrow(x) *
      (1 - (sum(diag(spline_est$hat)) / nrow(x)))^2

    smoothing_error[smoothing_idx] <- RSS / gcv_denom

    if (smoothing_error[smoothing_idx] < best_error) {
      best_error <- smoothing_error[smoothing_idx]
      best_spline <- spline_est
      opt_run <- smoothing_idx
    }

    # if errors have been increasing for past 4 smoothing values, stop early
    if (smoothing_idx >= 4) {
      recent_errors <- smoothing_error[(smoothing_idx - 3):smoothing_idx]
      if (all(recent_errors == cummax(recent_errors))) {
        break
      }
    }
  }

  embedding_map <- function(parameters) {
    as.vector(
      (t(best_spline$s) %*% etaFunc(parameters, params, 4 - d)) +
        (t(best_spline$alpha) %*% matrix(c(1, parameters), ncol = 1))
    )
  }

  return(
    list(
      embedding_map = embedding_map,
      hat = best_spline$hat,
      coefs = best_spline,
      lambda = smoothing_vals[opt_run]
    )
  )
}
