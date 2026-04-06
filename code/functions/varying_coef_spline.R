varying_coef_spline <- function(
  centers,
  params,
  weights,
  times,
  ids,
  scans,
  lambda,
  gamma,
  param_grid = NULL,
  verbose = FALSE
) {
  require(furrr, quietly = TRUE)
  require(here, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(Rfast, quietly = TRUE)

  source(here("code/functions/fit_weighted_spline.R"))

  d <- ncol(params) - 1
  D <- ncol(centers) - 1

  if (is.null(param_grid)) {
    # create new grid of parameter values, with number of points equal to maximum
    # of clusters across all scans
    N_prime <- max(table(scans))
    param_bounds <- colMinsMaxs(params[, -1])

    param_list <- list()
    for (dim in seq_len(d)) {
      param_range <- abs(param_bounds[2, dim] - param_bounds[1, dim])
      param_list[[dim]] <- seq(
        from = param_bounds[1, dim] - (0.1 * param_range),
        to = param_bounds[2, dim] + (0.1 * param_range),
        length.out = ceiling(N_prime^(1 / d))
      )
    }

    param_grid <- expand.grid(param_list)
  } else {
    param_grid <- as.matrix(param_grid)
  }

  grid_E <- calcE(param_grid, 4 - d)
  grid_mat <- cbind(rep(1, nrow(param_grid)), param_grid)

  n_knots <- nrow(param_grid)
  time_points <- sort(unique(times))

  n_scans <- vector()
  init_spline_est <- list()
  spline_est <- list()
  spline_coefs <- list()

  if (verbose == TRUE) {
    p <- progressor(along = time_points)
  }

  for (time_idx in seq_along(time_points)) {
    time_val <- time_points[time_idx]

    init_spline <- fit_weighted_spline(
      centers[times == time_val, -1],
      params[times == time_val, -1],
      weights[times == time_val],
      lambda
    )

    centers_grid <- map(
      seq_len(nrow(param_grid)),
      ~ init_spline$embedding_map(unlist(param_grid[.x, ]))
    ) |>
      reduce(rbind)

    spline_coefs[[time_idx]] <- solve_spline(
      grid_E,
      grid_mat,
      centers_grid,
      init_spline$lambda,
      d,
      D
    ) |>
      t() |>
      matrix(nrow = 1)

    if (verbose == TRUE) {
      p(message = sprintf("Time point %i of %i", time_idx, length(time_points)))
    }
  }

  spline_coefs <- reduce(spline_coefs, rbind)
  D_coef <- dim(spline_coefs)[2]

  time_mat <- cbind(rep(1, length(time_points)), time_points)
  time_E <- calcE(matrix(time_points, ncol = 1), 4 - 1)
  coef_gcv <- vector(mode = "numeric", length = length(gamma))

  best_time_spline <- NULL
  best_error <- 1e10
  opt_gamma <- 1

  for (smoothing_idx in seq_along(gamma)) {
    time_spline <- solve_spline_hat(
      time_E,
      time_mat,
      spline_coefs,
      gamma[smoothing_idx],
      1,
      D_coef
    )

    fitted_coefs <- time_spline$hat %*% spline_coefs

    coef_RSS <- apply(spline_coefs - fitted_coefs, 1, norm_euclidean)^2 |>
      sum()
    gcv_denom <- length(time_points) *
      (1 -
        (sum(diag(time_spline$hat)) / length(time_points)))^2

    coef_gcv[smoothing_idx] <- coef_RSS / gcv_denom

    if (coef_gcv[smoothing_idx] < best_error) {
      best_error <- coef_gcv[smoothing_idx]
      best_time_spline <- time_spline
      opt_run <- smoothing_idx
    }

    # if errors have been increasing for past 4 smoothing values, stop early
    if (smoothing_idx >= 4) {
      recent_errors <- coef_gcv[(smoothing_idx - 3):smoothing_idx]
      if (all(recent_errors == cummax(recent_errors))) {
        break
      }
    }
  }

  get_time_spline_coefs <- function(time_val) {
    return_vec <- as.vector(
      t(best_time_spline$s) %*%
        etaFunc(
          as.matrix(time_val),
          matrix(time_points, ncol = 1),
          3
        ) +
        t(best_time_spline$alpha) %*%
          matrix(c(1, time_val), ncol = 1)
    )
    return_vec
  }

  f_out <- function(param_vec) {
    coefs <- get_time_spline_coefs(param_vec[1])
    coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
    return_vec <- t(coef_mat[1:n_knots, ]) %*%
      etaFunc(param_vec[-1], param_grid, 4 - d) +
      t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*%
        matrix(c(1, param_vec[-1]), ncol = 1)
    c(param_vec[1], return_vec)
  }

  list(
    embedding_map = f_out,
    temporal_spline = best_time_spline,
    gamma = gamma[opt_gamma]
  )
}
