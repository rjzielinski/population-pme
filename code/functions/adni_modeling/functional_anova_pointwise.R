functional_anova_pointwise <- function(
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
  require(purrr, quietly = TRUE)

  n_partitions <- length(additive_model)
  n_groups <- length(additive_model[[1]]$group_embeddings)
  n_individuals <- length(additive_model[[1]]$id_embeddings)

  param_grids <- calc_param_grids(params, n_params, interval = 0.25)

  print("Computing Embeddings")

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

  test_stat <- list()
  rejected <- list()

  for (partition_idx in seq_len(n_partitions)) {
    test_stat[[partition_idx]] <- (sum_squares$ssh[[partition_idx]] /
      (n_groups - 1)) /
      (sum_squares$sse[[partition_idx]] / (n_individuals - n_groups))
  }

  if (bootstrap == FALSE) {
    if (test_type == "f_type") {
      critical_value <- qf(1 - alpha, n_groups - 1, n_individuals - n_groups)
    } else if (test_type == "chisq_type") {
      critical_value <- qchisq(1 - alpha, n_groups - 1) / (n_groups - 1)
    }
    for (partition_idx in seq_along(additive_model)) {
      rejected[[partition_idx]] <- test_stat[[partition_idx]] > critical_value
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

        test_stats <- map(
          seq_len(n_partitions),
          ~ (bootstrap_sum_squares$ssh[[.x]] /
            (n_groups - 1)) /
            (bootstrap_sum_squares$sse[[.x]] /
              (n_individuals - n_groups))
        )

        p()

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

calc_param_grids <- function(params, n_params, interval = 0.25) {
  require(Rfast, quietly = TRUE)

  param_bounds <- list()
  param_grids <- list()
  time_points_list <- list()

  for (param_idx in seq_along(params)) {
    time_points_list[[param_idx]] <- unique(params[[param_idx]][, 1])
    times <- seq(
      from = 0,
      to = max(time_points_list[[param_idx]]),
      by = interval
    )

    param_bounds[[param_idx]] <- colMinsMaxs(params[[param_idx]][, -1])
    d <- ncol(param_bounds[[param_idx]])

    param_list <- list()
    param_list[[1]] <- times

    for (dim in seq_len(d)) {
      param_range <- abs(
        param_bounds[[param_idx]][2, dim] - param_bounds[[param_idx]][1, dim]
      )
      param_list[[dim + 1]] <- seq(
        from = param_bounds[[param_idx]][2, dim] - (0.1 * param_range),
        to = param_bounds[[param_idx]][2, dim] + (0.1 * param_range),
        length.out = ceiling(n_params^(1 / d))
      )
    }

    param_grids[[param_idx]] <- expand.grid(param_list)
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
    p <- progressor(nrow(param_grids[[partition_idx]]))
    population_embeddings[[partition_idx]] <- future_map(
      seq_len(nrow(param_grids[[partition_idx]])),
      ~ {
        p()
        additive_model[[
          partition_idx
        ]]$population_embedding$embedding_map(unlist(param_grids[[
          partition_idx
        ]][
          .x,
        ]))
      },
      .options = furrr_options(seed = TRUE)
    ) |>
      reduce(rbind)

    part_group_embeddings <- foreach(
      group_idx = seq_len(n_groups),
      .options.future = list(seed = TRUE)
    ) %do%
      {
        p <- progressor(nrow(param_grids[[partition_idx]]))
        future_map(
          seq_len(nrow(param_grids[[partition_idx]])),
          ~ {
            p()
            additive_model[[partition_idx]]$group_embeddings[[
              group_idx
            ]]$embedding_map(unlist(param_grids[[partition_idx]][.x, ]))
          },
          .options = furrr_options(seed = TRUE)
        ) |>
          reduce(rbind)
      }

    p <- progressor(n_individuals)
    part_id_embeddings <- foreach(id_idx = seq_len(n_individuals)) %do%
      {
        p()
        future_map(
          seq_len(nrow(param_grids[[partition_idx]])),
          ~ additive_model[[partition_idx]]$id_embeddings[[
            id_idx
          ]]$embedding_map(unlist(param_grids[[partition_idx]][.x, ])),
          .options = furrr_options(seed = TRUE)
        ) |>
          reduce(rbind)
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

calc_sample_cov <- function() {
  # sample_cov <- list()
  # sample_cov[[partition_idx]] <- matrix(
  #   nrow = nrow(param_grids[[partition_idx]]),
  #   ncol = nrow(param_grids[[partition_idx]])
  # )
  #
  # for (param_idx1 in nrow(param_grids[[partition_idx]])) {
  #   for (param_idx2 in nrow(param_grids[[partition_idx]])) {
  #     sample_cov_num <- map(
  #       seq_along(ids),
  #       ~ sum(
  #         id_embeddings[[partition_idx]][[.x]][param_idx1, -1] *
  #           id_embeddings[[partition_idx]][[.x]][param_idx2, -1]
  #       )
  #     ) |>
  #       reduce(sum)
  #
  #     sample_cov[[partition_idx]][param_idx1, param_idx2] <- sample_cov_num /
  #       (n_individuals - n_groups)
  #   }
  # }
}
