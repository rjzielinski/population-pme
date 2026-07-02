functional_permutation_test <- function(
  additive_model,
  centers,
  weights,
  params,
  groups,
  ids,
  scans,
  times,
  partitions,
  id_groups,
  n_params,
  template,
  alpha = 0.05,
  n_permutations = 1000,
  threads = 4,
  mode = "additive_embeddings",
  progress = TRUE,
  verbose = TRUE
) {
  require(doFuture, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(future, quietly = TRUE)
  require(Matrix, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(RhpcBLASctl, quietly = TRUE)

  total_cores <- availableCores() - 4
  workers <- floor(total_cores / threads)

  old_plan <- plan(multicore, workers = workers)
  on.exit(plan(old_plan), add = TRUE)

  n_partitions <- length(additive_model)
  n_groups <- length(additive_model[[1]]$group_embeddings)
  n_individuals <- length(additive_model[[1]]$id_embeddings)

  group_vals <- sort(unique(groups))
  id_vals <- sort(unique(ids))
  scan_vals <- sort(unique(scans))
  partition_vals <- sort(unique(partitions))

  param_grids <- calc_param_grids(
    params,
    n_params,
    interval = 0.25,
    template = template
  )

  param_interval_vols <- vector()

  for (partition_idx in seq_len(n_partitions)) {
    param_interval_length <- vector()
    for (dim_idx in seq_len(ncol(param_grids[[partition_idx]]))) {
      values <- param_grids[[partition_idx]][, dim_idx] |>
        unique() |>
        sort()

      param_interval_length[dim_idx] <- values[2] - values[1]
    }

    param_interval_vols[partition_idx] <- prod(param_interval_length)
  }

  if (verbose == TRUE) {
    print("Computing Embeddings")
  }

  if (mode == "additive_embeddings") {
    embeddings <- calc_embeddings(
      additive_model,
      param_grids,
      n_groups,
      n_individuals,
      id_groups
    )

    sum_squares <- calc_sum_squares(
      embeddings,
      param_grids,
      groups,
      ids
    )
  } else if (mode == "sample_mean") {
    embeddings <- calc_embeddings_sample_means(
      additive_model,
      param_grids,
      n_groups,
      n_individuals,
      id_groups
    )

    sum_squares <- calc_sum_squares_sample_means(
      embeddings,
      param_grids,
      groups,
      ids
    )
  }

  gc()

  f_test_stat <- list()
  chisq_test_stat <- list()
  l2_norm_test_stat <- list()

  f_rejected <- list()
  chisq_rejected <- list()
  l2_norm_rejected <- list()

  f_p_values <- list()
  chisq_p_values <- list()
  l2_norm_p_values <- list()

  f_permute_rejected <- list()
  chisq_permute_rejected <- list()
  l2_norm_permute_rejected <- list()

  f_permute_p_values <- list()
  chisq_permute_p_values <- list()
  l2_norm_permute_p_values <- list()

  for (partition_idx in seq_len(n_partitions)) {
    f_test_stat[[partition_idx]] <- (sum_squares$ssh[[partition_idx]] /
      (n_groups - 1)) /
      (sum_squares$sse[[partition_idx]] / (n_individuals - n_groups))

    chisq_test_stat[[partition_idx]] <- (sum_squares$ssh[[partition_idx]] /
      (n_groups - 1)) /
      (sum_squares$sse[[partition_idx]] / (n_individuals - n_groups))

    l2_norm_test_stat[[partition_idx]] <- colSums(
      sum_squares$ssh[[partition_idx]] * param_interval_vols[partition_idx]
    )
  }

  f_critical_value <- qf(1 - alpha, n_groups - 1, n_individuals - n_groups)
  chisq_critical_value <- qchisq(1 - alpha, n_groups - 1) / (n_groups - 1)

  trace_sample_cov <- list()
  trace_sample_cov2 <- list()

  C <- list()
  D <- list()

  l2_norm_critical_values <- list()
  for (partition_idx in seq_len(n_partitions)) {
    trace_sample_cov[[partition_idx]] <- vector()
    trace_sample_cov2[[partition_idx]] <- vector()

    covariance_array <- calc_sample_cov_arr(
      embeddings,
      param_grids,
      partition_idx,
      n_individuals,
      n_groups
    )

    for (dim_idx in seq_len(dim(covariance_array)[3])) {
      covariance_mat <- covariance_array[,, dim_idx]

      trace_sample_cov[[partition_idx]][dim_idx] <- (diag(covariance_mat) *
        param_interval_vols[partition_idx]) |>
        sum()

      trace_sample_cov2[[partition_idx]][dim_idx] <- sum(covariance_mat^2) *
        param_interval_vols[partition_idx]^2
    }

    C[[partition_idx]] <- trace_sample_cov2[[partition_idx]] /
      trace_sample_cov[[partition_idx]]

    D[[partition_idx]] <- (trace_sample_cov[[partition_idx]]^2) /
      trace_sample_cov2[[partition_idx]]

    l2_norm_critical_values[[partition_idx]] <- C[[partition_idx]] *
      qchisq(1 - alpha, (n_groups - 1) * D[[partition_idx]])
  }

  for (partition_idx in seq_len(n_partitions)) {
    f_rejected[[partition_idx]] <- f_test_stat[[partition_idx]] >
      f_critical_value

    f_p_values[[partition_idx]] <- 1 -
      pf(
        f_test_stat[[partition_idx]],
        n_groups - 1,
        n_individuals - n_groups
      )

    chisq_rejected[[partition_idx]] <- chisq_test_stat[[partition_idx]] >
      chisq_critical_value

    chisq_p_values[[partition_idx]] <- 1 -
      pchisq(
        (n_groups - 1) * chisq_test_stat[[partition_idx]],
        n_groups - 1
      )

    l2_norm_rejected[[partition_idx]] <- l2_norm_test_stat[[partition_idx]] >
      l2_norm_critical_values[[partition_idx]]

    l2_norm_p_values[[partition_idx]] <- 1 -
      pchisq(
        l2_norm_test_stat[[partition_idx]] / C[[partition_idx]],
        (n_groups - 1) * D[[partition_idx]]
      )
  }

  f_permute_test_stats <- f_test_stat
  chisq_permute_test_stats <- chisq_test_stat
  l2_norm_permute_test_stats <- l2_norm_test_stat

  f_permute_test_stats <- map(
    f_permute_test_stats,
    ~ replicate(n_permutations, matrix(NA, nrow = nrow(.x), ncol = ncol(.x)))
  )
  chisq_permute_test_stats <- map(
    chisq_permute_test_stats,
    ~ replicate(n_permutations, matrix(NA, nrow = nrow(.x), ncol = ncol(.x)))
  )
  l2_norm_permute_test_stats <- map(
    l2_norm_permute_test_stats,
    ~ matrix(nrow = n_permutations, ncol = length(.x))
  )

  if (verbose == TRUE) {
    print("Getting Permutation Statistics...")
  }

  with_progress({
    if (progress == TRUE) {
      p <- progressor(n_permutations)
    }

    permutation_results <- foreach(
      permute_idx = seq_len(n_permutations),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        blas_set_num_threads(threads)
        omp_set_num_threads(threads)

        permute_id_groups <- id_groups |>
          sample(size = n_individuals)

        permute_groups <- vector(length = length(groups))
        for (id_idx in seq_along(id_vals)) {
          id_set <- ids == id_vals[id_idx]
          id_group_val <- permute_id_groups[id_idx]
          permute_groups[id_set] <- id_group_val
        }

        permute_model <- list()
        for (partition_idx in seq_len(n_partitions)) {
          partition_groups <- permute_groups[
            partitions == partition_vals[partition_idx]
          ]
          partition_ids <- ids[partitions == partition_vals[partition_idx]]
          partition_scans <- scans[partitions == partition_vals[partition_idx]]
          partition_times <- times[partitions == partition_vals[partition_idx]]

          permute_model[[partition_idx]] <- fit_permutation_models(
            additive_model[[partition_idx]]$population_embedding,
            centers[[partition_idx]],
            params[[partition_idx]],
            weights[[partition_idx]],
            additive_model[[partition_idx]]$param_grid,
            lambda = exp(-10:10),
            gamma = exp(-10:10),
            groups = partition_groups,
            ids = partition_ids,
            scans = partition_scans,
            times = partition_times,
            template = template,
            epsilon = 0.025,
            max_iter = 25,
            cores = 1,
            plot_progress = FALSE
          )
        }

        if (mode == "additive_embeddings") {
          permute_embeddings <- calc_embeddings(
            permute_model,
            param_grids,
            n_groups,
            n_individuals,
            permute_id_groups
          )

          gc()

          permute_sum_squares <- calc_sum_squares(
            permute_embeddings,
            param_grids,
            permute_groups,
            ids
          )
        } else if (mode == "sample_mean") {
          permute_embeddings <- calc_embeddings_sample_means(
            permute_model,
            param_grids,
            n_groups,
            n_individuals,
            permute_id_groups
          )

          gc()

          permute_sum_squares <- calc_sum_squares_sample_means(
            permute_embeddings,
            param_grids,
            permute_groups,
            ids
          )
        }

        permute_f_stat <- list()
        permute_chisq_stat <- list()
        permute_l2_norm_stat <- list()

        for (partition_idx in seq_len(n_partitions)) {
          permute_f_stat[[partition_idx]] <- (permute_sum_squares$ssh[[
            partition_idx
          ]] /
            (n_groups - 1)) /
            (permute_sum_squares$sse[[partition_idx]] /
              (n_individuals - n_groups))

          permute_chisq_stat[[partition_idx]] <- (permute_sum_squares$ssh[[
            partition_idx
          ]] /
            (n_groups - 1)) /
            (permute_sum_squares$sse[[partition_idx]] /
              (n_individuals - n_groups))

          permute_l2_norm_stat[[partition_idx]] <- colSums(
            permute_sum_squares$ssh[[partition_idx]] *
              param_interval_vols[partition_idx]
          )
        }

        if (progress == TRUE) {
          p(message = sprintf("Permutation sample: %g", permute_idx))
        }
        gc()

        list(
          f_test_stat = permute_f_stat,
          chisq_test_stat = permute_chisq_stat,
          l2_norm_test_stat = permute_l2_norm_stat
        )
      }
  })

  for (permute_idx in seq_len(n_permutations)) {
    for (partition_idx in seq_len(n_partitions)) {
      f_permute_test_stats[[partition_idx]][,,
        permute_idx
      ] <- permutation_results[[permute_idx]]$f_test_stat[[partition_idx]]

      chisq_permute_test_stats[[partition_idx]][,,
        permute_idx
      ] <- permutation_results[[permute_idx]]$chisq_test_stat[[partition_idx]]

      l2_norm_permute_test_stats[[partition_idx]][
        permute_idx,
      ] <- permutation_results[[permute_idx]]$l2_norm_test_stat[[partition_idx]]
    }
  }

  # calculate critical values at each parameter point
  f_permute_critical_values <- map(
    f_permute_test_stats,
    ~ matrix(
      data = NA,
      nrow = dim(.x)[1],
      ncol = dim(.x)[2]
    )
  )
  f_permute_p_vals <- map(
    f_permute_test_stats,
    ~ matrix(
      data = NA,
      nrow = dim(.x)[1],
      ncol = dim(.x)[2]
    )
  )

  chisq_permute_critical_values <- map(
    chisq_permute_test_stats,
    ~ matrix(
      data = NA,
      nrow = dim(.x)[1],
      ncol = dim(.x)[2]
    )
  )
  chisq_permute_p_vals <- map(
    chisq_permute_test_stats,
    ~ matrix(
      data = NA,
      nrow = dim(.x)[1],
      ncol = dim(.x)[2]
    )
  )

  l2_norm_permute_critical_values <- map(
    l2_norm_permute_test_stats,
    ~ vector(mode = "numeric", length = ncol(.x))
  )
  l2_norm_permute_p_vals <- map(
    l2_norm_permute_test_stats,
    ~ vector(mode = "numeric", length = ncol(.x))
  )

  for (partition_idx in seq_len(n_partitions)) {
    nrows <- dim(f_permute_test_stats[[partition_idx]])[1]
    ncols <- dim(f_permute_test_stats[[partition_idx]])[2]

    for (row_idx in seq_len(nrows)) {
      for (col_idx in seq_len(ncols)) {
        f_permute_critical_values[[partition_idx]][
          row_idx,
          col_idx
        ] <- quantile(
          f_permute_test_stats[[partition_idx]][row_idx, col_idx, ],
          probs = 1 - alpha,
          na.rm = TRUE
        )
        f_permute_p_vals[[partition_idx]][row_idx, col_idx] <- 1 -
          mean(
            f_test_stat[[partition_idx]][row_idx, col_idx] >
              f_permute_test_stats[[partition_idx]][row_idx, col_idx, ],
            na.rm = TRUE
          )

        chisq_permute_critical_values[[partition_idx]][
          row_idx,
          col_idx
        ] <- quantile(
          chisq_permute_test_stats[[partition_idx]][row_idx, col_idx, ],
          probs = 1 - alpha,
          na.rm = TRUE
        )
        chisq_permute_p_vals[[partition_idx]][row_idx, col_idx] <- 1 -
          mean(
            chisq_test_stat[[partition_idx]][row_idx, col_idx] >
              chisq_permute_test_stats[[partition_idx]][row_idx, col_idx, ],
            na.rm = TRUE
          )
      }
    }

    for (col_idx in seq_len(ncols)) {
      l2_norm_permute_critical_values[[partition_idx]][col_idx] <- quantile(
        l2_norm_permute_test_stats[[partition_idx]][, col_idx],
        probs = 1 - alpha,
        na.rm = TRUE
      )
      l2_norm_permute_p_vals[[partition_idx]][col_idx] <- 1 -
        mean(
          l2_norm_test_stat[[partition_idx]][col_idx] >
            l2_norm_permute_test_stats[[partition_idx]][, col_idx],
          na.rm = TRUE
        )
    }
  }

  for (partition_idx in seq_len(n_partitions)) {
    f_permute_rejected[[partition_idx]] <- f_test_stat[[partition_idx]] >
      f_permute_critical_values[[partition_idx]]
    f_permute_p_values[[partition_idx]] <- f_permute_p_vals[[partition_idx]]

    chisq_permute_rejected[[partition_idx]] <- chisq_test_stat[[
      partition_idx
    ]] >
      chisq_permute_critical_values[[partition_idx]]
    chisq_permute_p_values[[partition_idx]] <- chisq_permute_p_vals[[
      partition_idx
    ]]

    l2_norm_permute_rejected[[partition_idx]] <- l2_norm_test_stat[[
      partition_idx
    ]] >
      l2_norm_permute_critical_values[[partition_idx]]
    l2_norm_permute_p_values[[partition_idx]] <- l2_norm_permute_p_vals[
      partition_idx
    ]
  }

  f_test <- list(
    test_statistics = f_test_stat,
    param_critical_value = f_critical_value,
    param_rejected = f_rejected,
    param_p_values = f_p_values,
    permute_critical_values = f_permute_critical_values,
    permute_rejected = f_permute_rejected,
    permute_p_values = f_permute_p_values
  )

  chisq_test <- list(
    test_statistics = chisq_test_stat,
    param_critical_values = chisq_critical_value,
    param_rejected = chisq_rejected,
    param_p_values = chisq_p_values,
    permute_critical_values = chisq_permute_critical_values,
    permute_rejected = chisq_permute_rejected,
    permute_p_values = chisq_permute_p_values
  )

  l2_norm_test <- list(
    test_statistic = l2_norm_test_stat,
    param_critical_values = l2_norm_critical_values,
    param_rejected = l2_norm_rejected,
    param_p_value = l2_norm_p_values,
    permute_critical_values = l2_norm_permute_critical_values,
    permute_rejected = l2_norm_permute_rejected,
    permute_p_values = l2_norm_permute_p_values
  )

  test_out <- list(
    params = param_grids,
    embeddings = embeddings,
    l2_norm_test = l2_norm_test,
    f_test = f_test,
    chisq_test = chisq_test
  )
  test_out
}

calc_param_grids <- function(
  params,
  n_params,
  interval = 0.25,
  template = "euclidean",
  param_threshold = 0.0025
) {
  require(Rfast, quietly = TRUE)

  time_points_list <- list()
  param_bound_list <- list()
  param_grid_list <- list()
  param_interval_list <- list()

  for (param_idx in seq_along(params)) {
    time_points_list[[param_idx]] <- sort(unique(params[[param_idx]][, 1]))

    d <- ncol(params[[param_idx]]) - 1
    if (template == "euclidean") {
      param_bounds <- apply(
        params[[param_idx]][, -1, drop = FALSE],
        2,
        quantile,
        probs = c(param_threshold, 1 - param_threshold)
      )
      param_bound_list[[param_idx]] <- param_bounds
    }
  }

  all_times <- unlist(time_points_list) |>
    unique()
  min_time <- map(time_points_list, min) |> unlist() |> min()
  max_time <- map(time_points_list, max) |> unlist() |> max()
  n_times <- ceiling((max_time / interval) - (min_time / interval))
  time_vals <- seq(from = min_time, by = interval, length.out = n_times)

  param_grids <- list()
  for (param_idx in seq_along(params)) {
    param_grid_list <- list()
    if (template == "euclidean") {
      param_bound_mat <- param_bound_list[[param_idx]]
      param_intervals <- vector(mode = "numeric", length = d)

      param_list <- list()
      for (dim_idx in seq_len(d)) {
        param_intervals[dim_idx] <- (param_bound_mat[2, dim_idx] -
          param_bound_mat[1, dim_idx]) /
          n_params

        param_list[[dim]] <- seq(
          from = param_bound_mat[1, dim_idx],
          to = param_bound_mat[2, dim_idx],
          by = param_intervals[dim_idx]
        )
      }
      param_grid_list[[param_idx]] <- expand.grid(param_list) |>
        as.matrix()

      param_interval_list[[param_idx]] <- param_intervals
    } else if (template == "sphere") {
      param_grid_list[[param_idx]] <- fibonacci_sphere(n_params)
    }

    param_grids[[param_idx]] <- map(
      time_vals,
      ~ cbind(.x, param_grid_list[[param_idx]])
    ) |>
      reduce(rbind)
  }

  param_grids
}

calc_param_grids_resample <- function(
  params,
  n_params,
  interval = 0.25,
  noise_factor = 0,
  dist_kernel_sd = 1,
  template = "euclidean"
) {
  # manifolds of parameter values are not necessarily regularly shaped
  # manifolds may shift in shape and location over time
  # parameters are not uniformly distributed over the manifold
  #
  # instead of generating constant parameter grid, resample from
  # existing parameter values
  # sampling probabilities will be proportional to proximity between time points
  # add slight noise to sampled parameter values

  require(Rfast, quietly = TRUE)

  param_bounds <- list()
  param_grids <- list()
  time_points_list <- list()

  for (param_idx in seq_along(params)) {
    time_points_list[[param_idx]] <- sort(unique(params[[param_idx]][, 1]))

    time_params <- table(params[[param_idx]][, 1])

    time_idx_vals <- map(
      params[[param_idx]][, 1],
      ~ which(time_points_list[[param_idx]] == .x)
    ) |>
      reduce(c)

    times <- seq(
      from = 0,
      to = max(time_points_list[[param_idx]]),
      by = interval
    )

    param_list <- list()

    for (time_idx in seq_along(times)) {
      time_val <- times[time_idx]

      time_dists <- time_points_list[[param_idx]] - time_val

      dist_weights <- dnorm(time_dists, mean = 0, sd = dist_kernel_sd)
      dist_weights <- dist_weights / sum(dist_weights)

      time_weights <- dist_weights * (1 / time_params)

      param_weights <- time_weights[time_idx_vals]

      sample_idx <- sample(
        seq_len(nrow(params[[param_idx]])),
        size = n_params,
        prob = param_weights
      )

      sampled_params <- params[[param_idx]][sample_idx, ]
      sampled_params_jitter <- jitter(sampled_params, noise_factor)

      param_list[[time_idx]] <- cbind(time_val, sampled_params_jitter[, -1])
    }

    param_grids[[param_idx]] <- reduce(param_list, rbind)
  }

  param_grids
}

calc_embeddings <- function(
  additive_model,
  param_grids,
  n_groups,
  n_individuals,
  id_groups
) {
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(progressr, quietly = TRUE)

  population_embeddings <- list()
  group_embeddings <- list()
  id_embeddings <- list()

  grand_mean <- list()
  group_means <- list()

  for (partition_idx in seq_along(additive_model)) {
    param_grid <- param_grids[[partition_idx]]

    pop_embedding_map <- additive_model[[
      partition_idx
    ]]$population_embedding$embedding_map
    group_embedding_maps <- map(
      additive_model[[partition_idx]]$group_embeddings,
      ~ .x$embedding_map
    )
    id_embedding_maps <- map(
      additive_model[[partition_idx]]$id_embeddings,
      ~ .x$embedding_map
    )

    part_population_embeddings <- map(
      seq_len(nrow(param_grid)),
      ~ {
        pop_embedding_map(param_grid[.x, ])
      }
    )
    population_embeddings[[partition_idx]] <- do.call(
      rbind,
      part_population_embeddings
    )

    part_group_embeddings <- foreach(
      group_idx = seq_len(n_groups),
      .options.future = list(seed = TRUE)
    ) %do%
      {
        group_embedding_map <- group_embedding_maps[[group_idx]]
        group_embedding_list <- map(
          seq_len(nrow(param_grid)),
          ~ {
            group_embedding_map(param_grid[.x, ])
          }
        )
        do.call(rbind, group_embedding_list)
      }

    part_id_embeddings <- foreach(
      id_idx = seq_len(n_individuals),
      .options.future = list(seed = TRUE)
    ) %do%
      {
        id_embedding_map <- id_embedding_maps[[id_idx]]
        id_embedding_list <- map(
          seq_len(nrow(param_grid)),
          ~ id_embedding_map(param_grid[.x, ])
        )
        do.call(rbind, id_embedding_list)
      }

    full_embeddings <- map(
      seq_along(part_id_embeddings),
      ~ cbind(
        part_group_embeddings[[id_groups[.x]]][, 1],
        part_group_embeddings[[id_groups[.x]]][, -1] +
          part_id_embeddings[[.x]][, -1]
      )
    )

    grand_mean[[partition_idx]] <- reduce(full_embeddings, `+`) /
      length(full_embeddings)

    group_means[[partition_idx]] <- foreach(
      group_idx = seq_along(n_groups)
    ) %do%
      {
        full_embeddings[id_groups == group_idx] |>
          reduce(`+`) /
          sum(id_groups == group_idx)
      }

    group_embeddings[[partition_idx]] <- part_group_embeddings
    id_embeddings[[partition_idx]] <- part_id_embeddings
  }

  embedding_list <- list(
    population_embeddings = population_embeddings,
    group_embeddings = group_embeddings,
    id_embeddings = id_embeddings,
    group_assignments = id_groups
  )
}


calc_embeddings_sample_means <- function(
  additive_model,
  param_grids,
  n_groups,
  n_individuals,
  id_groups
) {
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(progressr, quietly = TRUE)

  full_embeddings <- list()
  population_embeddings <- list()
  group_embeddings <- list()
  id_embeddings <- list()

  grand_mean <- list()
  group_means <- list()

  for (partition_idx in seq_along(additive_model)) {
    param_grid <- param_grids[[partition_idx]]

    pop_embedding_map <- additive_model[[
      partition_idx
    ]]$population_embedding$embedding_map
    group_embedding_maps <- map(
      additive_model[[partition_idx]]$group_embeddings,
      ~ .x$embedding_map
    )
    id_embedding_maps <- map(
      additive_model[[partition_idx]]$id_embeddings,
      ~ .x$embedding_map
    )

    part_population_embeddings <- map(
      seq_len(nrow(param_grid)),
      ~ {
        pop_embedding_map(param_grid[.x, ])
      }
    )
    population_embeddings[[partition_idx]] <- do.call(
      rbind,
      part_population_embeddings
    )

    part_group_embeddings <- foreach(
      group_idx = seq_len(n_groups),
      .options.future = list(seed = TRUE)
    ) %do%
      {
        group_embedding_map <- group_embedding_maps[[group_idx]]
        group_embedding_list <- map(
          seq_len(nrow(param_grid)),
          ~ {
            group_embedding_map(param_grid[.x, ])
          }
        )
        do.call(rbind, group_embedding_list)
      }

    part_id_embeddings <- foreach(
      id_idx = seq_len(n_individuals),
      .options.future = list(seed = TRUE)
    ) %do%
      {
        id_embedding_map <- id_embedding_maps[[id_idx]]
        id_embedding_list <- map(
          seq_len(nrow(param_grid)),
          ~ id_embedding_map(param_grid[.x, ])
        )
        do.call(rbind, id_embedding_list)
      }

    full_embeddings[[partition_idx]] <- map(
      seq_along(part_id_embeddings),
      ~ cbind(
        part_group_embeddings[[id_groups[.x]]][, 1],
        population_embeddings[[partition_idx]][, -1] +
          part_group_embeddings[[id_groups[.x]]][, -1] +
          part_id_embeddings[[.x]][, -1]
      )
    )

    grand_mean[[partition_idx]] <- reduce(
      full_embeddings[[partition_idx]],
      `+`
    ) /
      length(full_embeddings[[partition_idx]])

    group_means[[partition_idx]] <- foreach(
      group_idx = seq_len(n_groups)
    ) %do%
      {
        full_embeddings[[partition_idx]][id_groups == group_idx] |>
          reduce(`+`) /
          sum(id_groups == group_idx)
      }

    # instead of using embeddings directly from embedding maps, use grand mean,
    # group-mean differences from grand mean, and differences between id embeddings
    # and group means

    group_embeddings[[partition_idx]] <- foreach(
      group_idx = seq_len(n_groups)
    ) %do%
      {
        cbind(
          grand_mean[[partition_idx]][, 1],
          group_means[[partition_idx]][[group_idx]][, -1] -
            grand_mean[[partition_idx]][, -1]
        )
      }

    id_embeddings[[partition_idx]] <- foreach(
      id_idx = seq_len(n_individuals)
    ) %do%
      {
        cbind(
          full_embeddings[[partition_idx]][[id_idx]][, 1],
          full_embeddings[[partition_idx]][[id_idx]][, -1] -
            group_means[[partition_idx]][[id_groups[id_idx]]][, -1]
        )
      }
  }

  embedding_list <- list(
    population_embeddings = population_embeddings,
    group_embeddings = group_embeddings,
    id_embeddings = id_embeddings,
    group_assignments = id_groups
  )
}


calc_sum_squares <- function(
  embeddings,
  param_grids,
  groups,
  ids
) {
  require(foreach, quietly = TRUE)

  population_embeddings <- embeddings$population_embeddings
  group_embeddings <- embeddings$group_embeddings
  id_embeddings <- embeddings$id_embeddings
  id_groups <- embeddings$group_assignments

  group_n <- map(groups, ~ sum(id_groups == .x)) |>
    reduce(c)

  ssh <- list()
  sse <- list()

  for (partition_idx in seq_along(param_grids)) {
    ssh[[partition_idx]] <- foreach(
      param_idx = seq_len(nrow(param_grids[[partition_idx]]))
    ) %do%
      {
        map(
          seq_along(group_embeddings[[partition_idx]]),
          ~ group_n[.x] *
            group_embeddings[[partition_idx]][[.x]][param_idx, -1]^2
        ) |>
          reduce(rbind) |>
          colSums()
      } |>
      reduce(rbind)

    sse[[partition_idx]] <- foreach(
      param_idx = seq_len(nrow(param_grids[[partition_idx]]))
    ) %do%
      {
        map(
          seq_along(id_embeddings[[partition_idx]]),
          ~ id_embeddings[[partition_idx]][[.x]][param_idx, -1]^2
        ) |>
          reduce(rbind) |>
          colSums()
      } |>
      reduce(rbind)
  }

  out <- list(
    ssh = ssh,
    sse = sse
  )
}

calc_sum_squares_sample_means <- function(
  embeddings,
  param_grids,
  groups,
  ids
) {
  require(foreach, quietly = TRUE)

  population_embeddings <- embeddings$population_embeddings
  group_embeddings <- embeddings$group_embeddings
  id_embeddings <- embeddings$id_embeddings
  id_groups <- embeddings$group_assignments

  group_n <- map(sort(unique(groups)), ~ sum(id_groups == .x)) |>
    reduce(c)

  n_partitions <- length(embeddings$population_embeddings)
  n_groups <- length(group_embeddings[[1]])
  n_individuals <- length(id_embeddings[[1]])

  covariance_arrays <- list()
  ssh <- list()
  sse <- list()

  for (partition_idx in seq_len(n_partitions)) {
    covariance_arrays[[partition_idx]] <- calc_sample_cov_arr(
      embeddings,
      param_grids,
      partition_idx,
      n_individuals,
      n_groups
    )

    ssh[[partition_idx]] <- map(
      seq_along(group_embeddings[[partition_idx]]),
      ~ group_n[.x] * group_embeddings[[partition_idx]][[.x]][, -1]^2
    ) |>
      reduce(`+`)

    sse[[partition_idx]] <- map(
      seq_len(dim(covariance_arrays[[partition_idx]])[3]),
      ~ (n_individuals - n_groups) *
        diag(covariance_arrays[[partition_idx]][,, .x])
    ) |>
      reduce(cbind)
  }

  out <- list(
    ssh = ssh,
    sse = sse
  )
}

calc_sample_cov_arr <- function(
  embeddings,
  param_grids,
  partition,
  n_individuals,
  n_groups
) {
  id_embeddings <- embeddings$id_embeddings[[partition]]
  D <- ncol(id_embeddings[[1]]) - 1

  complete_embeddings <- lapply(
    id_embeddings,
    function(x) x[, -1, drop = FALSE]
  ) |>
    unlist()

  embedding_array <- array(
    complete_embeddings,
    dim = c(nrow(param_grids[[partition]]), D, n_individuals)
  )

  covariance_array <- array(
    dim = c(dim(embedding_array)[1], dim(embedding_array)[1], D)
  )
  for (dim_idx in seq_len(D)) {
    dim_embedding <- embedding_array[, dim_idx, ]
    dim_cov_mat <- dim_embedding %*%
      t(dim_embedding) /
      (n_individuals - n_groups)
    covariance_array[,, dim_idx] <- dim_cov_mat
  }

  covariance_array
}

calc_sample_cov_mat <- function(
  embeddings,
  partition,
  n_individuals,
  n_groups
) {
  id_embeddings <- embeddings$id_embeddings[[partition]]

  # obtain embedding matrix for all IDs
  complete_embeddings <- lapply(
    id_embeddings,
    function(x) x[, -1, drop = FALSE]
  ) |>
    reduce(cbind)

  covariance_mat <- complete_embeddings %*%
    t(complete_embeddings) /
    (n_individuals - n_groups)

  covariance_mat
}

calc_sample_cov <- function(
  embeddings,
  param_grids,
  partition,
  param_idx1,
  param_idx2
) {
  # NOTE: this function assumes that param1 and param2 are elements of the param_grids matrices

  n_groups <- length(embeddings$group_embeddings[[partition]])
  id_embeddings <- embeddings$id_embeddings[[partition]]
  n_individuals <- length(id_embeddings)

  embeddings_param1 <- sapply(
    id_embeddings,
    function(x) x[param_idx1, -1]
  )

  if (param_idx1 == param_idx2) {
    sample_cov <- sum(embeddings_param1^2)
  } else {
    embeddings_param2 <- sapply(
      id_embeddings,
      function(x) x[param_idx2, -1]
    )

    sample_cov <- sum(embeddings_param1 * embeddings_param2)
  }

  sample_cov <- sample_cov / (n_individuals - n_groups)

  sample_cov
}

calc_sample_cov_vec <- function(
  embeddings,
  param_grids,
  partition,
  param_idx1,
  param_idx2
) {
  n_groups <- length(embeddings$group_embeddings[[partition]])
  id_embeddings <- embeddings$id_embeddings[[partition]]
  n_individuals <- length(id_embeddings)

  embeddings_param1 <- sapply(
    id_embeddings,
    function(x) x[param_idx1, -1]
  ) |>
    t()

  if (param_idx1 == param_idx2) {
    sample_cov <- colSums(embeddings_param1^2)
  } else {
    embeddings_param2 <- sapply(
      id_embeddings,
      function(x) x[param_idx2, -1]
    ) |>
      t()

    sample_cov <- colSums(embeddings_param1 * embeddings_param2)
  }

  sample_cov <- sample_cov / (n_individuals - n_groups)
  sample_cov
}

calc_sample_cov2 <- function(
  embeddings,
  param_grids,
  partition,
  param_idx1,
  param_idx2,
  param_interval_vols
) {
  sample_covs1 <- map_dbl(
    seq_len(nrow(param_grids[[partition]])),
    ~ calc_sample_cov(
      embeddings,
      param_grids,
      partition,
      param_idx1,
      .x
    )
  )

  if (param_idx1 == param_idx2) {
    sample_cov2_val <- sum(sample_covs1^2 * param_interval_vols[partition])
  } else {
    sample_covs2 <- map_dbl(
      seq_len(nrow(param_grids[[partition]])),
      ~ calc_sample_cov(
        embeddings,
        param_grids,
        partition,
        .x,
        param_idx2
      )
    )
    sample_cov2_val <- sum(
      sample_covs1 * sample_covs2 * param_interval_vols[partition]
    )
  }

  sample_cov2_val
}

calc_sample_cov2_opt <- function(
  embeddings,
  param_grids,
  partition,
  param_idx1,
  param_idx2,
  param_interval_vols
) {
  id_embeddings <- embeddings$id_embeddings[[partition]]
  n_individuals <- length(id_embeddings)
  n_groups <- length(embeddings$group_embeddings[[partition]])
  vol <- param_interval_vols[partition]

  # obtain embedding matrix for all IDs
  complete_embeddings <- lapply(
    id_embeddings,
    function(x) x[, -1, drop = FALSE]
  ) |>
    reduce(cbind)

  # tcrossprod() is a faster version of mat %*% t(mat)
  covariance_mat <- tcrossprod(complete_embeddings) / (n_individuals - n_groups)

  sample_covs1 <- covariance_mat[param_idx1, ]

  if (param_idx1 == param_idx2) {
    sample_cov2_val <- sum(sample_covs1^2 * vol)
  } else {
    sample_covs2 <- covariance_mat[param_idx2, ]
    sample_cov2_val <- sum(sample_covs1 * sample_covs2 * vol)
  }

  sample_cov2_val
}


fit_permutation_models <- function(
  population_embedding,
  centers,
  params,
  weights,
  param_grid,
  lambda,
  gamma,
  groups,
  ids,
  scans,
  times,
  template,
  epsilon = 0.05,
  max_iter = 100,
  cores = 1,
  verbose = FALSE,
  plot_progress = FALSE
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
  source(here("code/functions/gen_full_embedding.R"))
  source(here("code/functions/varying_coef_spline.R"))

  D <- ncol(centers) - 1
  d <- ncol(params) - 1

  group_vals <- sort(unique(groups))
  id_vals <- sort(unique(ids))
  scan_vals <- sort(unique(scans))

  if (verbose == TRUE) {
    print("Initializing...")
  }

  # create new grid of parameter values, with number of points equal to maximum
  # of clusters across all scans
  if (is.null(param_grid)) {
    N_prime <- max(table(scans))

    if (template == "euclidean") {
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
    } else if (template == "sphere") {
      param_grid <- fibonacci_sphere(N_prime)
    }
  }

  population_preds <- map(
    seq_len(nrow(params)),
    ~ population_embedding$embedding_map(params[.x, ])
  ) %>%
    reduce(rbind)

  group_preds <- matrix(0, nrow = nrow(centers), ncol = ncol(centers))
  id_preds <- matrix(0, nrow = nrow(centers), ncol = ncol(centers))

  if (verbose == TRUE) {
    print("Estimating initial group-level embeddings...")
  }

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
        template,
        param_grid = param_grid,
        verbose = verbose
      )
    }
  names(group_embeddings) <- group_vals

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

  if (verbose == TRUE) {
    print("Estimating initial individual-level embeddings...")
  }

  id_embeddings <- foreach(
    id_idx = seq_along(id_vals),
    .options.future = list(seed = TRUE)
  ) %do%
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

      varying_coef_spline(
        adjusted_centers,
        id_params,
        id_weights,
        id_times,
        rep(id_vals[id_idx], nrow(id_centers)),
        id_scans,
        lambda,
        gamma,
        template,
        param_grid = param_grid
      )
    }
  names(id_embeddings) <- id_vals

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

  if (plot_progress == TRUE) {
    if (D == 3) {
      fig <- plot_ly(
        x = full_preds[, 2],
        y = full_preds[, 3],
        z = full_preds[, 4],
        frame = full_preds[, 1],
        color = groups,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)
      )
      print(fig)
    } else if (D == 2) {
      fig <- plot_ly(
        x = full_preds[, 2],
        y = full_preds[, 1],
        z = full_preds[, 3],
        color = groups,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)
      )
      print(fig)
    }
  }

  mse <- map(
    seq_len(nrow(centers)),
    ~ weights[.x] * dist_euclidean(centers[.x, ], full_preds[.x, ])^2
  ) %>%
    reduce(c) %>%
    mean()

  if (verbose == TRUE) {
    print(
      paste0(
        "Backfitting Iteration ",
        as.character(0),
        ": ",
        "Estimated mean squared error - ",
        as.character(round(mse, 10)),
        "; Relative change in mean squared error - ",
        as.character(NA)
      )
    )
  }

  if (verbose == TRUE) {
    print("Initialization Complete, beginning iterations")
  }

  while ((epsilon_hat > epsilon) & (n <= max_iter)) {
    mse_old <- mse
    if (verbose == TRUE) {
      print("Estimating group-level embeddings...")
    }

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
          template,
          param_grid = param_grid,
          verbose = verbose
        )
      }
    names(group_embeddings) <- group_vals

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

    if (verbose == TRUE) {
      print("Estimating individual-level embeddings...")
    }

    id_embeddings <- foreach(
      id_idx = seq_along(id_vals),
      .options.future = list(seed = TRUE)
    ) %do%
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

        varying_coef_spline(
          adjusted_centers,
          id_params,
          id_weights,
          id_times,
          rep(id_vals[id_idx], nrow(id_centers)),
          id_scans,
          lambda,
          gamma,
          template,
          param_grid = param_grid
        )
      }
    names(id_embeddings) <- id_vals

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

    if (plot_progress == TRUE) {
      if (D == 3) {
        fig <- plot_ly(
          x = full_preds[, 2],
          y = full_preds[, 3],
          z = full_preds[, 4],
          frame = full_preds[, 1],
          color = groups,
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 3)
        )
        print(fig)
      } else if (D == 2) {
        fig <- plot_ly(
          x = full_preds[, 2],
          y = full_preds[, 1],
          z = full_preds[, 3],
          color = groups,
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 3)
        )
        print(fig)
      }
    }

    mse <- map(
      seq_len(nrow(centers)),
      ~ weights[.x] * dist_euclidean(centers[.x, ], full_preds[.x, ])^2
    ) %>%
      reduce(c) %>%
      mean()

    mse_ratio <- abs(mse - mse_old) / mse_old
    epsilon_hat <- mse_ratio
    n <- n + 1

    if (verbose == TRUE) {
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
  }

  full_embeddings <- foreach(id_idx = seq_along(id_vals)) %do%
    {
      id_set <- ids == id_vals[id_idx]
      group_idx <- which(group_vals == unique(groups[id_set]))

      gen_full_embedding(
        population_embedding,
        group_embeddings,
        id_embeddings,
        group_idx,
        id_idx
      )
    }

  return(
    list(
      embeddings = full_embeddings,
      population_embedding = population_embedding,
      group_embeddings = group_embeddings,
      id_embeddings = id_embeddings,
      param_grid = param_grid,
      group_values = group_vals,
      id_values = id_vals
    )
  )
}
