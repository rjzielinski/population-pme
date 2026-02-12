fit_adni_additive_model <- function(
  centers,
  params,
  weights,
  lambda,
  gamma,
  k,
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
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(Rfast, quietly = TRUE)

  source(here("code/functions/fit_weighted_spline.R"))
  source(here("code/functions/varying_coef_spline.R"))

  avail_cores <- cores %/% k

  D <- ncol(centers) - 1
  d <- ncol(params) - 1

  group_vals <- unique(groups)[order(unique(groups))]
  id_vals <- unique(ids)[order(unique(ids))]
  scan_vals <- unique(scans)[order(unique(scans))]

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

  init_population_embedding <- varying_coef_spline(
    centers,
    params,
    weights,
    times,
    ids,
    scans,
    lambda,
    gamma,
    k,
    param_grid = param_grid,
    verbose = TRUE
  )

  print("Estimating initial group-level embeddings...")

  init_group_embeddings <- list()
  group_centers <- list()
  group_params <- list()
  group_weights <- list()

  # should I use dofuture and set seed?
  group_out <- foreach(group_idx = seq_along(group_vals)) %do%
    {
      group_set <- groups == group_vals[group_idx]
      temp_group_centers <- centers[group_set, ]
      temp_group_params <- params[group_set, ]
      temp_group_weights <- weights[group_set]
      temp_group_ids <- ids[group_set]
      temp_group_times <- times[group_set]
      temp_group_scans <- scans[group_set]

      if (length(unique(temp_group_ids)) < k) {
        temp_k <- length(unique(temp_group_ids))
      } else {
        temp_k <- k
      }

      temp_group_preds <- map(
        seq_len(nrow(temp_group_params)),
        ~ init_population_embedding$embedding_map(temp_group_params[.x, ])
      ) %>%
        reduce(rbind)

      adjusted_centers <- cbind(
        temp_group_centers[, 1],
        temp_group_centers[, -1] - temp_group_preds[, -1]
      )

      init_group_embedding <- varying_coef_spline(
        adjusted_centers,
        temp_group_params,
        temp_group_weights,
        temp_group_times,
        temp_group_ids,
        temp_group_scans,
        lambda,
        gamma,
        k,
        param_grid = param_grid,
        verbose = TRUE
      )

      list(
        embedding = init_group_embedding,
        group_centers = temp_group_centers,
        group_params = temp_group_params,
        group_weights = temp_group_weights,
        group_times = temp_group_times,
        group_ids = temp_group_ids,
        group_scans = temp_group_scans
      )
    }

  for (group_idx in seq_along(group_vals)) {
    init_group_embeddings[[group_idx]] <- group_out[[group_idx]]$embedding
    group_centers[[group_idx]] <- group_out[[group_idx]]$group_centers
    group_params[[group_idx]] <- group_out[[group_idx]]$group_params
    group_weights[[group_idx]] <- group_out[[group_idx]]$group_weights
  }

  init_id_embeddings <- list()
  id_centers <- list()
  id_params <- list()
  id_weights <- list()

  print("Estimating initial individual-level embeddings...")

  p <- progressor(along = unique(ids))
  id_out <- foreach(
    id_idx = seq_along(id_vals),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      id_set <- ids == id_vals[id_idx]
      temp_id_centers <- centers[id_set, ]
      temp_id_params <- params[id_set, ]
      temp_id_weights <- weights[id_set]
      temp_id_times <- times[id_set]
      temp_id_scans <- scans[id_set]
      temp_id_group <- unique(groups[id_set])
      group_idx <- which(group_vals == temp_id_group)

      temp_id_preds <- map(
        seq_len(nrow(temp_id_params)),
        ~ init_population_embedding$embedding_map(temp_id_params[.x, ]) +
          init_group_embeddings[[group_idx]]$embedding_map(temp_id_params[.x, ])
      ) %>%
        reduce(rbind)
      temp_id_preds <- cbind(temp_id_params[, 1], temp_id_preds[, -1])

      adjusted_centers <- cbind(
        temp_id_centers[, 1],
        temp_id_centers[, -1] - temp_id_preds[, -1]
      )

      init_id_embedding <- varying_coef_spline(
        adjusted_centers,
        temp_id_params,
        temp_id_weights,
        temp_id_times,
        rep(id_vals[id_idx], nrow(temp_id_centers)),
        temp_id_scans,
        lambda,
        gamma,
        k,
        param_grid = param_grid
      )

      p(message = sprintf("ID: %s", id_vals[id_idx]))
      list(
        embedding = init_id_embedding,
        id_centers = temp_id_centers,
        id_params = temp_id_params,
        id_weights = temp_id_weights,
        id_times = temp_id_times,
        id_group = temp_id_group,
        id_scans = temp_id_scans
      )
    }

  for (id_idx in seq_along(id_vals)) {
    init_id_embeddings[[id_idx]] <- id_out[[id_idx]]$embedding
    id_centers[[id_idx]] <- id_out[[id_idx]]$id_centers
    id_params[[id_idx]] <- id_out[[id_idx]]$id_params
    id_weights[[id_idx]] <- id_out[[id_idx]]$id_weights
  }

  epsilon_hat <- 2 * epsilon
  n <- 0
  population_embedding <- init_population_embedding
  group_embeddings <- init_group_embeddings
  id_embeddings <- init_id_embeddings

  center_preds <- map(
    seq_len(nrow(centers)),
    ~ population_embedding$embedding_map(params[.x, ]) +
      group_embeddings[[which(group_vals == groups[.x])]]$embedding_map(params[
        .x,
      ]) +
      id_embeddings[[which(id_vals == ids[.x])]]$embedding_map(params[.x, ])
  ) %>%
    reduce(rbind)
  center_preds <- cbind(centers[, 1], center_preds[, -1])

  mse <- map(
    seq_len(nrow(centers)),
    ~ weights[.x] * dist_euclidean(centers[.x, ], center_preds[.x, ])^2
  ) %>%
    reduce(c) %>%
    mean()

  print("Initialization Complete, beginning iterations")

  while ((epsilon_hat > epsilon) & (n <= max_iter)) {
    population_embedding_old <- population_embedding
    group_embeddings_old <- group_embeddings
    id_embeddings_old <- id_embeddings
    mse_old <- mse

    centers_pred_nopop <- map(
      seq_len(nrow(centers)),
      ~ group_embeddings[[which(
        group_vals == groups[.x]
      )]]$embedding_map(params[.x, ]) +
        id_embeddings[[which(id_vals == ids[.x])]]$embedding_map(params[
          .x,
        ])
    ) %>%
      reduce(rbind)
    centers_pred_nopop <- cbind(centers[, 1], centers_pred_nopop[, -1])

    adjusted_centers <- cbind(
      centers[, 1],
      centers[, -1] - centers_pred_nopop[, -1]
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
      k,
      param_grid = param_grid,
      verbose = TRUE
    )

    print("Estimating group-level embeddings...")

    group_out <- foreach(group_idx = seq_along(group_vals)) %do%
      {
        group_set <- groups == group_vals[group_idx]
        temp_group_centers <- centers[group_set, ]
        temp_group_params <- params[group_set, ]
        temp_group_weights <- weights[group_set]
        temp_group_ids <- ids[group_set]
        temp_group_times <- times[group_set]
        temp_group_scans <- scans[group_set]

        if (length(unique(temp_group_ids)) < k) {
          temp_k <- length(unique(temp_group_ids))
        } else {
          temp_k <- k
        }

        temp_group_preds <- map(
          seq_len(nrow(temp_group_params)),
          ~ (population_embedding$embedding_map(temp_group_params[.x, ]) +
            id_embeddings[[which(
              id_vals == temp_group_ids[.x]
            )]]$embedding_map(temp_group_params[.x, ]))
        ) %>%
          reduce(rbind)

        adjusted_centers <- cbind(
          temp_group_centers[, 1],
          temp_group_centers[, -1] - temp_group_preds[, -1]
        )

        group_embedding <- varying_coef_spline(
          adjusted_centers,
          temp_group_params,
          temp_group_weights,
          temp_group_times,
          temp_group_ids,
          temp_group_scans,
          lambda,
          gamma,
          k,
          param_grid = param_grid,
          verbose = TRUE
        )

        list(
          embedding = group_embedding,
          group_centers = temp_group_centers,
          group_params = temp_group_params,
          group_weights = temp_group_weights,
          group_times = temp_group_times,
          group_ids = temp_group_ids,
          group_scans = temp_group_scans
        )
      }

    for (group_idx in seq_along(unique(groups))) {
      group_embeddings[[group_idx]] <- group_out[[group_idx]]$embedding
      group_centers[[group_idx]] <- group_out[[group_idx]]$group_centers
      group_params[[group_idx]] <- group_out[[group_idx]]$group_params
      group_weights[[group_idx]] <- group_out[[group_idx]]$group_weights
    }

    print("Estimating individual-level embeddings...")

    p <- progressor(along = id_vals)
    id_out <- foreach(
      id_idx = seq_along(id_vals),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        id_set <- ids == id_vals[id_idx]
        temp_id_centers <- centers[id_set, ]
        temp_id_params <- params[id_set, ]
        temp_id_weights <- weights[id_set]
        temp_id_times <- times[id_set]
        temp_id_scans <- scans[id_set]
        temp_id_group <- unique(groups[id_set])
        group_idx <- which(unique(groups) == temp_id_group)

        temp_id_preds <- map(
          seq_len(nrow(temp_id_params)),
          ~ population_embedding$embedding_map(temp_id_params[.x, ]) +
            group_embeddings[[
              group_idx
            ]]$embedding_map(temp_id_params[.x, ])
        ) %>%
          reduce(rbind)

        adjusted_centers <- cbind(
          temp_id_centers[, 1],
          temp_id_centers[, -1] - temp_id_preds[, -1]
        )

        id_embedding <- varying_coef_spline(
          adjusted_centers,
          temp_id_params,
          temp_id_weights,
          temp_id_times,
          rep(id_vals[id_idx], nrow(temp_id_centers)),
          temp_id_scans,
          lambda,
          gamma,
          k,
          param_grid = param_grid
        )

        p(message = sprintf("ID: %s", id_vals[id_idx]))

        list(
          embedding = id_embedding,
          id_centers = temp_id_centers,
          id_params = temp_id_params,
          id_weights = temp_id_weights,
          id_times = temp_id_times,
          id_group = temp_id_group,
          id_scans = temp_id_scans
        )
      }

    for (id_idx in seq_along(unique(ids))) {
      id_embeddings[[id_idx]] <- id_out[[id_idx]]$embedding
      id_centers[[id_idx]] <- id_out[[id_idx]]$id_centers
      id_params[[id_idx]] <- id_out[[id_idx]]$id_params
      id_weights[[id_idx]] <- id_out[[id_idx]]$id_weights
    }

    n <- n + 1

    center_preds <- map(
      seq_len(nrow(centers)),
      ~ (population_embedding$embedding_map(params[.x, ]) +
        group_embeddings[[which(
          group_vals == groups[.x]
        )]]$embedding_map(params[.x, ]) +
        id_embeddings[[which(id_vals == ids[.x])]]$embedding_map(params[
          .x,
        ]))
    ) %>%
      reduce(rbind)
    center_preds <- cbind(centers[, 1], center_preds[, -1])

    mse <- map(
      seq_len(nrow(centers)),
      ~ weights[.x] * dist_euclidean(centers[.x, ], center_preds[.x, ])^2
    ) %>%
      reduce(c) %>%
      mean()

    mse_ratio <- abs(mse - mse_old) / mse_old

    epsilon_hat <- mse_ratio

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
