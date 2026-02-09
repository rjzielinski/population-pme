varying_coef_spline <- function(
  centers,
  params,
  weights,
  times,
  ids,
  scans,
  lambda,
  gamma,
  k,
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
  }

  param_grid <- as.matrix(param_grid)
  n_knots <- nrow(param_grid)
  time_points <- unique(times)
  time_points <- time_points[order(time_points)]

  n_scans <- vector()
  init_spline_est <- list()
  spline_est <- list()
  spline_coefs <- list()

  if (verbose == TRUE) {
    p <- progressor(along = time_points)
  }

  for (time_idx in seq_along(time_points)) {
    time_val <- time_points[time_idx]

    scan_vals <- scans[times == time_val]
    n_scans[time_idx] <- scan_vals |>
      unique() |>
      length()

    fold_vec_full <- rep(1:k, ceiling(n_scans[time_idx] / k))

    scan_folds <- fold_vec_full[
      sample(
        seq_len(n_scans[time_idx]),
        n_scans[time_idx],
        replace = FALSE
      )
    ]

    folds <- future_map(
      scan_vals,
      ~ scan_folds[which(unique(scan_vals) == .x)]
    ) |>
      reduce(c)

    init_spline_est[[time_idx]] <- fit_weighted_spline(
      centers[times == time_val, -1],
      params[times == time_val, -1],
      weights[times == time_val],
      lambda,
      folds
    )

    centers_grid <- map(
      seq_len(nrow(param_grid)),
      ~ init_spline_est[[time_idx]]$embedding_map(unlist(param_grid[.x, ]))
    ) |>
      reduce(rbind)

    grid_E <- calcE(param_grid, 4 - d)
    grid_mat <- cbind(rep(1, nrow(param_grid)), param_grid)

    grid_coefs <- solve_spline(
      grid_E,
      grid_mat,
      centers_grid,
      init_spline_est[[time_idx]]$lambda,
      d,
      D
    ) |>
      t() |>
      matrix(nrow = 1)

    spline_coefs[[time_idx]] <- grid_coefs

    embedding_map <- function(parameters) {
      as.vector(
        (t(spline_coefs[seq_len(nrow(centers_grid)), ]) %*%
          etaFunc(parameters, param_grid, 4 - d)) +
          (t(spline_coefs[
            (nrow(centers_grid) + 1):(nrow(centers_grid) + d + 1),
          ]) %*%
            matrix(c(1, parameters), ncol = 1))
      )
    }

    if (verbose == TRUE) {
      p(message = sprintf("Time point %i of %i", time_idx, length(time_points)))
    }

    spline_est[[time_idx]] <- list(
      embedding_map = embedding_map,
      coefs = grid_coefs,
      lambda = init_spline_est[[time_idx]]$lambda
    )
  }

  spline_coefs <- reduce(spline_coefs, rbind)
  D_coef <- dim(spline_coefs)[2]

  mse <- vector()

  for (smoothing_idx in seq_along(gamma)) {
    cv_folds <- sample(
      seq_along(time_points),
      length(time_points),
      replace = FALSE
    )

    cv_mse <- vector()
    for (fold_idx in seq_along(cv_folds)) {
      time_E <- calcE(matrix(time_points[-cv_folds[fold_idx]], ncol = 1), 4 - 1)
      time_mat <- cbind(
        rep(1, length(time_points) - 1),
        time_points[-cv_folds[fold_idx]]
      )

      time_coefs <- solve_spline(
        time_E,
        time_mat,
        spline_coefs[-cv_folds[fold_idx], ],
        gamma[smoothing_idx],
        1,
        D_coef
      )

      get_time_spline_coefs <- function(time_val) {
        return_vec <- as.vector(
          t(time_coefs[1:(nrow(spline_coefs) - 1), ]) %*%
            etaFunc(
              as.matrix(time_val),
              matrix(time_points[-cv_folds[fold_idx]], ncol = 1),
              3
            ) +
            t(time_coefs[
              (length(time_points)):(length(time_points) + 1),
            ]) %*%
              matrix(c(1, time_val), ncol = 1)
        )
        return_vec
      }

      f_new <- function(param_vec) {
        coefs <- get_time_spline_coefs(param_vec[1])
        coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
        return_vec <- t(coef_mat[1:n_knots, ]) %*%
          etaFunc(param_vec[-1], param_grid, 4 - d) +
          t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*%
            matrix(c(1, param_vec[-1]), ncol = 1)
        c(param_vec[1], return_vec)
      }

      fold_centers <- centers[centers[, 1] == time_points[cv_folds[fold_idx]], ]
      fold_params <- params[params[, 1] == time_points[cv_folds[fold_idx]], ]

      fold_preds <- map(
        seq_len(nrow(fold_params)),
        ~ f_new(unlist(fold_params[.x, ]))
      ) |>
        reduce(rbind)

      cv_mse[fold_idx] <- map(
        seq_len(nrow(fold_centers)),
        ~ dist_euclidean(
          unlist(fold_centers[.x, -1]),
          unlist(fold_preds[.x, -1])
        )^2
      ) |>
        reduce(c) |>
        mean()
    }
    mse[smoothing_idx] <- mean(cv_mse)
  }

  opt_gamma <- which.min(mse)

  time_E <- calcE(matrix(time_points, ncol = 1), 4 - 1)
  time_mat <- cbind(
    rep(1, length(time_points)),
    time_points
  )

  time_coefs <- solve_spline(
    time_E,
    time_mat,
    spline_coefs,
    gamma[opt_gamma],
    1,
    D_coef
  )

  get_time_spline_coefs <- function(time_val) {
    return_vec <- as.vector(
      t(time_coefs[seq_len(nrow(spline_coefs)), ]) %*%
        etaFunc(
          as.matrix(time_val),
          matrix(time_points, ncol = 1),
          3
        ) +
        t(time_coefs[
          (length(time_points) + 1):(length(time_points) + 1 + 1),
        ]) %*%
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
    coefs = time_coefs,
    gamma = gamma[opt_gamma]
  )
}
