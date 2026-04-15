functional_anova <- function(
  additive_model,
  params,
  groups,
  ids,
  id_groups,
  n_params,
  test_type = "f_type",
  alpha = 0.05,
  bootstrap = FALSE,
  n_bootstrap = 1000
) {
  require(doFuture, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(future, quietly = TRUE)
  require(Matrix, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  n_partitions <- length(additive_model)
  n_groups <- length(additive_model[[1]]$group_embeddings)
  n_individuals <- length(additive_model[[1]]$id_embeddings)

  param_grids <- calc_param_grids(
    params,
    n_params,
    interval = 0.25,
    noise_factor = 10,
    dist_kernel_sd = 0.5
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

  print("Computing Embeddings")

  embeddings <- calc_embeddings(
    additive_model,
    param_grids,
    n_groups,
    n_individuals,
    id_groups
  )

  gc()

  sum_squares <- calc_sum_squares(
    embeddings,
    param_grids,
    groups,
    ids
  )

  test_stat <- list()
  rejected <- list()

  for (partition_idx in seq_len(n_partitions)) {
    if (test_type %in% c("f_type", "chisq_type")) {
      test_stat[[partition_idx]] <- (sum_squares$ssh[[partition_idx]] /
        (n_groups - 1)) /
        (sum_squares$sse[[partition_idx]] / (n_individuals - n_groups))
    } else if (test_type == "l2_norm") {
      param_interval_length <- vector()

      test_stat[[partition_idx]] <- sum(
        sum_squares$ssh[[partition_idx]] * param_interval_vols[partition_idx]
      )
    }
  }

  if (bootstrap == FALSE) {
    if (test_type == "f_type") {
      critical_value <- qf(1 - alpha, n_groups - 1, n_individuals - n_groups)
    } else if (test_type == "chisq_type") {
      critical_value <- qchisq(1 - alpha, n_groups - 1) / (n_groups - 1)
    } else if (test_type == "l2_norm") {
      sample_cov <- list()
      trace_sample_cov <- list()
      trace_sample_cov2 <- list()
      C <- list()
      D <- list()

      critical_value <- list()
      for (partition_idx in seq_len(n_partitions)) {
        p <- progressor(nrow(param_grids[[partition_idx]]))
        trace_sample_cov[[partition_idx]] <- future_map(
          seq_len(nrow(param_grids[[partition_idx]])),
          ~ {
            p()
            calc_sample_cov(
              embeddings,
              param_grids,
              partition_idx,
              .x,
              .x
            ) *
              param_interval_vols[partition_idx]
          },
          .options = furrr_options(seed = TRUE)
        ) |>
          reduce(sum)

        p <- progressor(nrow(param_grids[[partition_idx]]))
        trace_sample_cov2[[partition_idx]] <- future_map(
          seq_len(nrow(param_grids[[partition_idx]])),
          ~ {
            p()
            calc_sample_cov2(
              embeddings,
              param_grids,
              partition_idx,
              .x,
              .x,
              param_interval_vols
            ) *
              param_interval_vols[partition_idx]
          },
          .options = furrr_options(seed = TRUE)
        ) |>
          reduce(sum)

        C[[partition_idx]] <- trace_sample_cov2[[partition_idx]] /
          trace_sample_cov[[partition_idx]]

        D[[partition_idx]] <- (trace_sample_cov[[partition_idx]]^2) /
          trace_sample_cov2[[partition_idx]]

        critical_value[[partition_idx]] <- C[[partition_idx]] *
          qchisq(1 - alpha, (n_groups - 1) * D[[partition_idx]])
      }
    }

    for (partition_idx in seq_along(additive_model)) {
      if (test_type == "l2_norm") {
        rejected[[partition_idx]] <- test_stat[[partition_idx]] >
          critical_value[[partition_idx]]
      } else {
        rejected[[partition_idx]] <- test_stat[[partition_idx]] > critical_value
      }
    }
  } else {
    bootstrap_test_stats <- test_stat

    bootstrap_test_stats <- map(
      bootstrap_test_stats,
      ~ matrix(NA, nrow = n_bootstrap, ncol = length(.x))
    )

    print("Getting Bootstrap Sample...")

    p <- progressor(n_bootstrap)

    bootstrap_results <- foreach(
      bootstrap_idx = seq_len(n_bootstrap),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        bootstrap_embeddings <- calc_bootstrap_embeddings(
          embeddings,
          groups,
          ids,
          id_groups
        )

        bootstrap_sum_squares <- calc_sum_squares(
          bootstrap_embeddings,
          param_grids,
          groups,
          ids
        )

        if (test_type %in% c("f_type", "chisq_type")) {
          test_stats <- map(
            seq_len(n_partitions),
            ~ (bootstrap_sum_squares$ssh[[.x]] / (n_groups - 1)) /
              (bootstrap_sum_squares$sse[[.x]] / (n_individuals - n_groups))
          )
        } else if (test_type == "l2_norm") {
          test_stats <- map(
            seq_len(n_partitions),
            ~ sum(
              bootstrap_sum_squares$ssh[[partition_idx]] *
                param_interval_vols[.x]
            )
          )
        }

        p(message = sprintf("Bootstrap sample: %g", bootstrap_idx))
        gc()

        test_stats
      }

    for (bootstrap_idx in seq_len(n_bootstrap)) {
      for (partition_idx in seq_len(n_partitions)) {
        bootstrap_test_stats[[partition_idx]][
          bootstrap_idx,
        ] <- bootstrap_results[[bootstrap_idx]][[partition_idx]]
      }
    }

    # calculate critical values at each parameter point
    critical_value <- map(
      bootstrap_test_stats,
      ~ apply(.x, 2, quantile, 1 - alpha, na.rm = TRUE)
    )

    for (partition_idx in seq_len(n_partitions)) {
      rejected[[partition_idx]] <- test_stat[[partition_idx]] >
        critical_value[[partition_idx]]
    }
  }

  test_out <- list(
    params = param_grids,
    embeddings = embeddings,
    rejected = rejected
  )
}

calc_param_grids <- function(
  params,
  n_params,
  interval = 0.25,
  noise_factor = 0,
  dist_kernel_sd = 1
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

    with_progress({
      p <- progressor(n_groups)
      part_group_embeddings <- foreach(
        group_idx = seq_len(n_groups),
        .options.future = list(seed = TRUE)
      ) %dofuture%
        {
          group_embedding_map <- group_embedding_maps[[group_idx]]
          p(sprintf("Group %d of %d", group_idx, n_groups))
          group_embedding_list <- map(
            seq_len(nrow(param_grid)),
            ~ {
              group_embedding_map(param_grid[.x, ])
            }
          )
          do.call(rbind, group_embedding_list)
        }
    })

    with_progress({
      p <- progressor(n_individuals)
      part_id_embeddings <- foreach(
        id_idx = seq_len(n_individuals),
        .options.future = list(seed = TRUE)
      ) %dofuture%
        {
          id_embedding_map <- id_embedding_maps[[id_idx]]
          p(sprintf("ID %d of %d", id_idx, n_individuals))
          id_embedding_list <- map(
            seq_len(nrow(param_grid)),
            ~ id_embedding_map(param_grid[.x, ])
          )
          do.call(rbind, id_embedding_list)
        }
    })

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
            sum(group_embeddings[[partition_idx]][[.x]][param_idx, -1]^2)
        ) |>
          reduce(sum)
      } |>
      reduce(c)

    sse[[partition_idx]] <- foreach(
      param_idx = seq_len(nrow(param_grids[[partition_idx]]))
    ) %do%
      {
        map(
          seq_along(id_embeddings[[partition_idx]]),
          ~ sum(id_embeddings[[partition_idx]][[.x]][param_idx, -1]^2)
        ) |>
          reduce(sum)
      } |>
      reduce(c)
  }

  out <- list(
    ssh = ssh,
    sse = sse
  )
}

calc_bootstrap_embeddings <- function(embeddings, groups, ids, id_groups) {
  population_embeddings <- embeddings$population_embeddings
  group_embeddings <- embeddings$group_embeddings
  id_embeddings <- embeddings$id_embeddings

  bootstrap_id_idx <- sample(seq_along(ids), size = length(ids), replace = TRUE)
  bootstrap_groups <- id_groups[bootstrap_id_idx] |>
    sample(size = length(bootstrap_id_idx), replace = FALSE)

  bootstrap_group_n <- map(groups, ~ sum(bootstrap_groups == .x)) |>
    reduce(c)

  bootstrap_embeddings <- list()
  bootstrap_pop_embeddings <- list()
  bootstrap_group_embeddings <- list()
  bootstrap_id_embeddings <- list()

  for (partition_idx in seq_along(population_embeddings)) {
    partition_embeddings <- map(
      bootstrap_id_idx,
      ~ {
        time_vals <- population_embeddings[[partition_idx]][, 1]
        embedding_vals <- cbind(
          time_vals,
          population_embeddings[[partition_idx]][, -1] +
            group_embeddings[[partition_idx]][[which(
              groups == id_groups[.x]
            )]][, -1] +
            id_embeddings[[partition_idx]][[.x]][, -1]
        )
        embedding_vals
      }
    )

    bootstrap_embeddings[[partition_idx]] <- partition_embeddings

    bootstrap_pop_embeddings[[partition_idx]] <- reduce(
      bootstrap_embeddings[[partition_idx]],
      `+`
    ) /
      length(bootstrap_embeddings[[partition_idx]])

    partition_group_embeddings <- foreach(group_idx = seq_along(groups)) %do%
      {
        group_indices <- which(bootstrap_groups == groups[group_idx])

        group_mean <- (reduce(
          bootstrap_embeddings[[partition_idx]][group_indices],
          `+`
        ) /
          bootstrap_group_n[group_idx])

        group_mean_adj <- group_mean - bootstrap_pop_embeddings[[partition_idx]]
      }
    bootstrap_group_embeddings[[partition_idx]] <- partition_group_embeddings

    partition_id_embeddings <- foreach(
      id_idx = seq_along(bootstrap_id_idx)
    ) %do%
      {
        group_idx <- which(groups == bootstrap_groups[id_idx])
        id_adj <- bootstrap_embeddings[[partition_idx]][[id_idx]] -
          bootstrap_pop_embeddings[[partition_idx]] -
          bootstrap_group_embeddings[[partition_idx]][[group_idx]]
      }
    bootstrap_id_embeddings[[partition_idx]] <- partition_id_embeddings
  }

  bootstrap_embeddings <- list(
    population_embeddings = bootstrap_pop_embeddings,
    group_embeddings = bootstrap_group_embeddings,
    id_embeddings = bootstrap_id_embeddings,
    group_assignments = bootstrap_groups
  )
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
  n_individuals <- length(embeddings$id_embeddings[[partition]])

  embeddings_param1 <- map(
    embeddings$id_embeddings[[partition]],
    ~ .x[param_idx1, ][-1]
  )
  embeddings_param1 <- do.call(rbind, embeddings_param1)

  if (param_idx1 == param_idx2) {
    sample_cov <- diag(embeddings_param1 %*% t(embeddings_param1)) |>
      sum()
  } else {
    embeddings_param2 <- map(
      embeddings$id_embeddings[[partition]],
      ~ .x[param_idx2, ][-1]
    )
    embeddings_param2 <- do.call(rbind, embeddings_param2)

    sample_cov <- diag(embeddings_param1 %*% t(embeddings_param2)) |>
      sum()
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
  sample_covs1 <- map(
    seq_len(nrow(param_grids[[partition]])),
    ~ calc_sample_cov(
      embeddings,
      param_grids,
      partition,
      param_idx1,
      .x
    )
  ) |>
    reduce(c)

  if (param_idx1 == param_idx2) {
    sample_cov2_val <- sum(sample_covs1^2 * param_interval_vols[partition])
  } else {
    sample_covs2 <- map(
      seq_len(nrow(param_grids[[partition]])),
      ~ calc_sample_cov(
        embeddings,
        param_grids,
        partition,
        .x,
        param_idx2
      )
    ) |>
      reduce(c)

    sample_cov2_val <- sum(
      sample_covs1 * sample_covs2 * param_interval_vols[partition]
    )
  }

  sample_cov2_val
}
