fit_weighted_spline <- function(x, params, weights, smoothing_vals) {
  require(furrr, quietly = TRUE)
  require(future, quietly = TRUE)
  require(magrittr, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(purrr, quietly = TRUE)

  coefs <- list()
  d <- ncol(params)
  smoothing_error <- vector()
  lambda <- smoothing_vals

  for (smoothing_idx in seq_along(smoothing_vals)) {
    train_param_mat <- cbind(rep(1, nrow(params)), params)
    train_weights <- diag(weights)
    train_E <- calcE(params, 4 - d)

    coefs[[smoothing_idx]] <- solve_weighted_spline_hat(
      train_E,
      train_weights,
      train_param_mat,
      x,
      lambda[smoothing_idx],
      ncol(params),
      ncol(x)
    )

    fitted_values <- coefs[[smoothing_idx]]$hat %*% x

    RSS <- apply(x - fitted_values, 1, norm_euclidean)^2 |>
      sum()

    gcv_denom <- nrow(x) *
      (1 - (sum(diag(coefs[[smoothing_idx]]$hat)) / nrow(x)))^2

    smoothing_error[smoothing_idx] <- RSS / gcv_denom
    #
    # if errors have been increasing for past 4 smoothing values, stop early
    if (smoothing_idx >= 4) {
      recent_errors <- smoothing_error[(smoothing_idx - 3):smoothing_idx]
      if (all(recent_errors == cummax(recent_errors))) {
        break
      }
    }
  }

  opt_run <- which.min(smoothing_error)

  opt_coefs <- coefs[[opt_run]]

  embedding_map <- function(parameters) {
    as.vector(
      (t(opt_coefs$s) %*% etaFunc(parameters, params, 4 - d)) +
        (t(opt_coefs$alpha) %*% matrix(c(1, parameters), ncol = 1))
    )
  }

  return(
    list(
      embedding_map = embedding_map,
      hat = opt_coefs$hat,
      coefs = opt_coefs,
      lambda = smoothing_vals[opt_run]
    )
  )
}
