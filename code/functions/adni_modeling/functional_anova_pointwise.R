functional_anova_pointwise <- function(
  additive_model,
  params,
  groups,
  ids,
  id_groups,
  n_params,
  test_type = "f_type",
  alpha = 0.05
) {
  require(furrr)
  require(purrr)

  param_bounds <- list()
  param_grids <- list()
  time_points_list <- list()

  n_partitions <- length(additive_model)
  n_groups <- length(additive_model[[1]]$group_embeddings)
  n_individuals <- length(additive_model[[1]]$id_embeddings)

  group_n <- map(groups, ~ sum(id_groups == .x)) |>
    reduce(c)

  population_embeddings <- list()
  group_embeddings <- list()
  id_embeddings <- list()

  for (param_idx in seq_len(length(params))) {
    time_points_list[[param_idx]] <- unique(params[[param_idx]][, 1])
    times <- seq(
      from = 0,
      to = max(time_points_list[[param_idx]]),
      by = 0.25
    )

    param_bounds[[param_idx]] <- colMinsMaxs(lhipp_test_out$params[[
      param_idx
    ]][,
      -1
    ])
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

  for (partition_idx in seq_along(additive_model)) {
    population_embeddings[[partition_idx]] <- future_map(
      seq_len(nrow(param_grids[[partition_idx]])),
      ~ additive_model[[
        partition_idx
      ]]$population_embedding$embedding_map(unlist(param_grids[[partition_idx]][
        .x,
      ])),
      .options = furrr_options(seed = TRUE)
    ) |>
      reduce(rbind)

    # with_progress({
    #   p <- progressor(n_groups)
    part_group_embeddings <- foreach(
      group_idx = seq_len(n_groups),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        map(
          seq_len(nrow(param_grids[[partition_idx]])),
          ~ additive_model[[partition_idx]]$group_embeddings[[
            group_idx
          ]]$embedding_map(unlist(param_grids[[partition_idx]][.x, ]))
        ) |>
          reduce(rbind)
      }
    # })

    # with_progress({
    #   p <- progressor(n_individuals)
    part_id_embeddings <- foreach(
      id_idx = seq_len(n_individuals),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        map(
          seq_len(nrow(param_grids[[partition_idx]])),
          ~ additive_model[[partition_idx]]$id_embeddings[[
            id_idx
          ]]$embedding_map(unlist(param_grids[[partition_idx]][.x, ]))
        ) |>
          reduce(rbind)
        #      p()
      }
    # })

    group_embeddings[[partition_idx]] <- part_group_embeddings
    id_embeddings[[partition_idx]] <- part_id_embeddings
  }

  ssh <- list()
  sse <- list()
  # sample_cov <- list()
  f_score <- list()
  chisq_score <- list()

  rejected <- list()
  p_vals <- list()

  for (partition_idx in seq_along(additive_model)) {
    ssh[[partition_idx]] <- foreach(
      param_idx = seq_len(nrow(param_grids[[partition_idx]]))
    ) %do%
      {
        map(
          seq_along(groups),
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
          seq_along(ids),
          ~ sum(id_embeddings[[partition_idx]][[.x]][param_idx, -1]^2)
        ) |>
          reduce(sum)
      } |>
      reduce(c)

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

    f_score[[partition_idx]] <- (ssh[[partition_idx]] / (n_groups - 1)) /
      (sse[[partition_idx]] / (n_individuals - n_groups))
  }

  if (test_type == "f_type") {
    critical_value <- qf(1 - alpha, n_groups - 1, n_individuals - n_groups)
  } else if (test_type == "chisq_type") {
    critical_value <- qchisq(1 - alpha, n_groups - 1) / (n_groups - 1)
  }

  for (partition_idx in seq_along(additive_model)) {
    rejected[[partition_idx]] <- f_score[[partition_idx]] > critical_value
  }

  test_out <- list(
    params = param_grids,
    rejected = rejected
  )
}
