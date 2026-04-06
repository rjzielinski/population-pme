fit_adni_additive_model <- function(
  centers,
  params,
  weights,
  lambda,
  gamma,
  groups,
  ids,
  scans,
  times,
  epsilon = 0.05,
  max_iter = 100,
  cores = 1
) {
  require(doFuture, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(future, quietly = TRUE)
  require(here, quietly = TRUE)
  require(plotly, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(Rfast, quietly = TRUE)

  source(here("code/functions/fit_weighted_spline.R"))
  source(here("code/functions/varying_coef_spline.R"))

  D <- ncol(centers) - 1
  d <- ncol(params) - 1

  group_vals <- sort(unique(groups))
  id_vals <- sort(unique(ids))
  scan_vals <- sort(unique(scans))

  print("Initializing...")

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

  print("Estimating initial population-level embedding...")

  population_embedding <- varying_coef_spline(
    centers,
    params,
    weights,
    times,
    ids,
    scans,
    lambda,
    gamma,
    param_grid = param_grid,
    verbose = TRUE
  )

  population_preds <- map(
    seq_len(nrow(params)),
    ~ population_embedding$embedding_map(params[.x, ])
  ) %>%
    reduce(rbind)

  group_preds <- matrix(0, nrow = nrow(centers), ncol = ncol(centers))
  id_preds <- matrix(0, nrow = nrow(centers), ncol = ncol(centers))

  print("Estimating initial group-level embeddings...")

  group_embeddings <- foreach(group_idx = seq_along(group_vals)) %do%
    {
      group_set <- groups == group_vals[group_idx]

      group_centers <- centers[group_set, ]
      group_params <- params[group_set, ]
      group_ids <- ids[group_set]
      group_weights <- weights[group_set]
      group_times <- times[group_set]
      group_scans <- scans[group_set]

      adjusted_centers <- cbind(
        group_centers[, 1],
        group_centers[, -1] - population_preds[group_set, -1]
      )

      varying_coef_spline(
        adjusted_centers,
        group_params,
        group_weights,
        group_times,
        group_ids,
        group_scans,
        lambda,
        gamma,
        param_grid = param_grid,
        verbose = TRUE
      )
    }

  for (group_idx in seq_along(group_vals)) {
    group_set <- which(groups == group_vals[group_idx])
    temp_group_params <- params[group_set, ]

    group_preds[group_set, ] <- apply(
      temp_group_params,
      1,
      group_embeddings[[group_idx]]$embedding_map
    ) |>
      t()
  }

  print("Estimating initial individual-level embeddings...")

  p <- progressor(along = unique(ids))

  id_embeddings <- foreach(
    id_idx = seq_along(id_vals),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      id_set <- ids == id_vals[id_idx]

      id_centers <- centers[id_set, ]
      id_params <- params[id_set, ]
      id_weights <- weights[id_set]
      id_times <- times[id_set]
      id_scans <- scans[id_set]
      id_group <- unique(groups[id_set])

      group_idx <- which(group_vals == id_group)

      adjusted_centers <- cbind(
        id_centers[, 1],
        id_centers[, -1] -
          (population_preds[id_set, -1] + group_preds[id_set, -1])
      )

      p(message = sprintf("ID: %s", id_vals[id_idx]))

      varying_coef_spline(
        adjusted_centers,
        id_params,
        id_weights,
        id_times,
        rep(id_vals[id_idx], nrow(id_centers)),
        id_scans,
        lambda,
        gamma,
        param_grid = param_grid
      )
    }

  for (id_idx in seq_along(id_vals)) {
    id_set <- which(ids == id_vals[id_idx])
    temp_id_params <- params[id_set, ]

    id_preds[id_set, ] <- apply(
      temp_id_params,
      1,
      id_embeddings[[id_idx]]$embedding_map
    ) |>
      t()
  }

  epsilon_hat <- 2 * epsilon
  n <- 0

  full_preds <- cbind(
    centers[, 1],
    population_preds[, -1] + group_preds[, -1] + id_preds[, -1]
  )

  mse <- map(
    seq_len(nrow(centers)),
    ~ weights[.x] * dist_euclidean(centers[.x, ], full_preds[.x, ])^2
  ) %>%
    reduce(c) %>%
    mean()

  print("Initialization Complete, beginning iterations")

  while ((epsilon_hat > epsilon) & (n <= max_iter)) {
    mse_old <- mse

    adjusted_centers <- cbind(
      centers[, 1],
      centers[, -1] - (group_preds[, -1] + id_preds[, -1])
    )

    print("Estimating population-level embedding...")

    population_embedding <- varying_coef_spline(
      adjusted_centers,
      params,
      weights,
      times,
      ids,
      scans,
      lambda,
      gamma,
      param_grid = param_grid,
      verbose = TRUE
    )

    print("Estimating group-level embeddings...")

    group_embeddings <- foreach(group_idx = seq_along(group_vals)) %do%
      {
        group_set <- groups == group_vals[group_idx]

        group_centers <- centers[group_set, ]
        group_params <- params[group_set, ]
        group_ids <- ids[group_set]
        group_weights <- weights[group_set]
        group_times <- times[group_set]
        group_scans <- scans[group_set]

        temp_preds <- population_preds[group_set, ]
        adjusted_centers <- cbind(
          group_centers[, 1],
          group_centers[, -1] -
            (population_preds[group_set, -1] + id_preds[group_set, -1])
        )

        varying_coef_spline(
          adjusted_centers,
          group_params,
          group_weights,
          group_times,
          group_ids,
          group_scans,
          lambda,
          gamma,
          param_grid = param_grid,
          verbose = TRUE
        )
      }

    for (group_idx in seq_along(group_vals)) {
      group_set <- which(groups == group_vals[group_idx])
      temp_group_params <- params[group_set, ]

      group_preds[group_set, ] <- apply(
        temp_group_params,
        1,
        group_embeddings[[group_idx]]$embedding_map
      ) |>
        t()
    }

    print("Estimating individual-level embeddings...")

    p <- progressor(along = unique(ids))

    id_embeddings <- foreach(
      id_idx = seq_along(id_vals),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        id_set <- ids == id_vals[id_idx]

        id_centers <- centers[id_set, ]
        id_params <- params[id_set, ]
        id_weights <- weights[id_set]
        id_times <- times[id_set]
        id_scans <- scans[id_set]
        id_group <- unique(groups[id_set])

        group_idx <- which(group_vals == id_group)

        adjusted_centers <- cbind(
          id_centers[, 1],
          id_centers[, -1] -
            (population_preds[id_set, -1] + group_preds[id_set, -1])
        )

        p(message = sprintf("ID: %s", id_vals[id_idx]))

        varying_coef_spline(
          adjusted_centers,
          id_params,
          id_weights,
          id_times,
          rep(id_vals[id_idx], nrow(id_centers)),
          id_scans,
          lambda,
          gamma,
          param_grid = param_grid
        )
      }

    for (id_idx in seq_along(id_vals)) {
      id_set <- which(ids == id_vals[id_idx])
      temp_id_params <- params[id_set, ]

      id_preds[id_set, ] <- apply(
        temp_id_params,
        1,
        id_embeddings[[id_idx]]$embedding_map
      ) |>
        t()
    }

    full_preds <- cbind(
      centers[, 1],
      population_preds[, -1] + group_preds[, -1] + id_preds[, -1]
    )

    mse <- map(
      seq_len(nrow(centers)),
      ~ weights[.x] * dist_euclidean(centers[.x, ], full_preds[.x, ])^2
    ) %>%
      reduce(c) %>%
      mean()

    mse_ratio <- abs(mse - mse_old) / mse_old
    epsilon_hat <- mse_ratio
    n <- n + 1

    print(
      paste0(
        "Backfitting Iteration ",
        as.character(n),
        ": ",
        "Estimated mean squared error - ",
        as.character(round(mse, 10)),
        "; Relative change in mean squared error - ",
        as.character(round(mse_ratio, 5))
      )
    )
  }

  full_embeddings <- list()
  for (id_idx in seq_along(id_vals)) {
    id_set <- ids == id_vals[id_idx]
    group_idx <- which(group_vals == unique(groups[id_set]))

    full_embeddings[[id_idx]] <- function(parameters) {
      pred_vec <- population_embedding$embedding_map(parameters) +
        group_embeddings[[group_idx]]$embedding_map(parameters) +
        id_embeddings[[id_idx]]$embedding_map(parameters)
      c(parameters[1], pred_vec[-1])
    }
  }

  return(
    list(
      embeddings = full_embeddings,
      population_embedding = population_embedding,
      group_embeddings = group_embeddings,
      id_embeddings = id_embeddings
    )
  )
}
