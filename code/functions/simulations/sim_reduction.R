sim_reduction <- function(
  sim_preprocessed,
  case,
  min_clusters,
  component_type,
  subsample_size = 5
) {
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)

  if (case == 1) {
    red_names <- c(
      "X1",
      "X2",
      "weight",
      "group",
      "id",
      "scan",
      "time_from_bl"
    )
  } else if (case == 10) {
    red_names <- c(
      "X1",
      "X2",
      "X3",
      "weight",
      "group",
      "id",
      "scan",
      "time_from_bl"
    )
  }

  population_scans <- sim_preprocessed$data_full |>
    filter(id == "Population") |>
    pull(scan) |>
    unique()

  group_scans <- foreach(group_val = sim_preprocessed$groups) %do%
    {
      sim_preprocessed$data_full |>
        filter(id == paste0("Group ", group_val)) |>
        pull(scan) |>
        unique()
    }

  with_progress({
    p <- progressor(steps = length(sim_preprocessed$scans))

    reduced_dfs <- foreach(
      scan_idx = seq_along(sim_preprocessed$scans),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        scan_df <- sim_preprocessed$data |>
          filter(scan == sim_preprocessed$scans[scan_idx])

        group_val <- unique(scan_df$group)
        id_val <- unique(scan_df$id)

        scan_data <- scan_df |>
          select(contains("X"), -contains("true")) |>
          as.matrix()

        scan_red <- hdmde(
          scan_data,
          min_clusters,
          0.05,
          200
        )

        scan_centers <- scan_red$mu
        scan_weights <- scan_red$theta_hat
        scan_n <- nrow(scan_centers)

        scan_time <- unique(scan_df$time)

        if (component_type == "centers") {
          scan_red <- cbind(
            scan_centers,
            scan_weights,
            rep(group_val, scan_n),
            rep(id_val, scan_n),
            rep(sim_preprocessed$scans[scan_idx], scan_n),
            rep(scan_time, scan_n)
          ) |>
            as_tibble(.name_repair = c("minimal"))
        } else if (component_type == "subsample") {
          km <- scan_red$km

          cluster_points <- foreach(cluster_val = seq_len(scan_n)) %do%
            {
              n_obs <- km$size[cluster_val]
              cluster_obs <- scan_data[
                km$cluster == cluster_val,
                ,
                drop = FALSE
              ]

              n_points <- min(n_obs, subsample_size)

              if (n_obs >= subsample_size) {
                sampled_indices <- sample(n_obs, subsample_size)
                cluster_points <- cluster_obs[sampled_indices, ]
                point_weights <- rep(
                  scan_weights[cluster_val] / subsample_size,
                  subsample_size
                )
              } else {
                cluster_points <- cluster_obs
                point_weights <- rep(
                  scan_weights[cluster_val] / n_obs,
                  n_obs
                )
              }

              cbind(
                cluster_points,
                point_weights,
                rep(group_val, n_points),
                rep(id_val, n_points),
                rep(sim_preprocessed$scans[scan_idx], n_points),
                rep(scan_time, n_points)
              )
            }

          scan_red <- do.call(rbind, cluster_points) |>
            as_tibble(.name_repair = c("minimal"))
        }

        names(scan_red) <- red_names

        p(sprintf("Scan #%i of %i", scan_idx, length(sim_preprocessed$scans)))

        scan_red
      }
  })

  sim_reduced <- do.call(bind_rows, reduced_dfs)

  if (case == 1) {
    sim_reduced <- sim_reduced |>
      mutate(
        X1 = as.numeric(X1),
        X2 = as.numeric(X2),
        weight = as.numeric(weight),
        time_from_bl = as.numeric(time_from_bl),
        partition = 1
      )
  } else if (case == 10) {
    sim_reduced <- sim_reduced |>
      mutate(
        X1 = as.numeric(X1),
        X2 = as.numeric(X2),
        X3 = as.numeric(X3),
        weight = as.numeric(weight),
        time_from_bl = as.numeric(time_from_bl),
        partition = 1
      )
  }

  with_progress({
    p <- progressor(steps = length(population_scans))

    population_reduced_dfs <- foreach(
      scan_idx = seq_along(population_scans),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        scan_df <- sim_preprocessed$data_full |>
          filter(scan == population_scans[scan_idx])

        group_val <- unique(scan_df$group)
        id_val <- unique(scan_df$id)

        scan_data <- scan_df |>
          select(contains("X"), -contains("true")) |>
          as.matrix()

        scan_red <- hdmde(
          scan_data,
          min_clusters,
          0.05,
          200
        )

        scan_centers <- scan_red$mu
        scan_weights <- scan_red$theta_hat
        scan_n <- nrow(scan_centers)

        scan_time <- unique(scan_df$time)

        if (component_type == "centers") {
          scan_red <- cbind(
            scan_centers,
            scan_weights,
            rep(group_val, scan_n),
            rep(id_val, scan_n),
            rep(population_scans[scan_idx], scan_n),
            rep(scan_time, scan_n)
          ) |>
            as_tibble(.name_repair = c("minimal"))
        } else if (component_type == "subsample") {
          km <- scan_red$km

          cluster_points <- foreach(cluster_val = seq_len(scan_n)) %do%
            {
              n_obs <- km$size[cluster_val]
              cluster_obs <- scan_data[
                km$cluster == cluster_val,
                ,
                drop = FALSE
              ]

              n_points <- min(n_obs, subsample_size)

              if (n_obs >= subsample_size) {
                sampled_indices <- sample(n_obs, subsample_size)
                cluster_points <- cluster_obs[sampled_indices, ]
                point_weights <- rep(
                  scan_weights[cluster_val] / subsample_size,
                  subsample_size
                )
              } else {
                cluster_points <- cluster_obs
                point_weights <- rep(
                  scan_weights[cluster_val] / n_obs,
                  n_obs
                )
              }

              cbind(
                cluster_points,
                point_weights,
                rep(group_val, n_points),
                rep(id_val, n_points),
                rep(population_scans[scan_idx], n_points),
                rep(scan_time, n_points)
              )
            }

          scan_red <- do.call(rbind, cluster_points) |>
            as_tibble(.name_repair = c("minimal"))
        }

        names(scan_red) <- red_names

        p(sprintf("Scan #%i of %i", scan_idx, length(population_scans)))

        scan_red
      }
  })

  population_reduced <- do.call(bind_rows, population_reduced_dfs)

  if (case == 1) {
    population_reduced <- population_reduced |>
      mutate(
        X1 = as.numeric(X1),
        X2 = as.numeric(X2),
        weight = as.numeric(weight),
        time_from_bl = as.numeric(time_from_bl),
        partition = 1
      )
  } else if (case == 10) {
    population_reduced <- population_reduced |>
      mutate(
        X1 = as.numeric(X1),
        X2 = as.numeric(X2),
        X3 = as.numeric(X3),
        weight = as.numeric(weight),
        time_from_bl = as.numeric(time_from_bl),
        partition = 1
      )
  }

  groups_reduced <- list()
  for (group_idx in seq_along(sim_preprocessed$groups)) {
    with_progress({
      p <- progressor(steps = length(group_scans[[group_idx]]))

      group_reduced_dfs <- foreach(
        scan_idx = seq_along(group_scans[[group_idx]]),
        .options.future = list(seed = TRUE)
      ) %dofuture%
        {
          scan_df <- sim_preprocessed$data_full |>
            filter(scan == group_scans[[group_idx]][scan_idx])

          group_val <- unique(scan_df$group)
          id_val <- unique(scan_df$id)

          scan_data <- scan_df |>
            select(contains("X"), -contains("true")) |>
            as.matrix()

          scan_red <- hdmde(
            scan_data,
            min_clusters,
            0.05,
            200
          )

          scan_centers <- scan_red$mu
          scan_weights <- scan_red$theta_hat
          scan_n <- nrow(scan_centers)

          scan_time <- unique(scan_df$time)

          if (component_type == "centers") {
            scan_red <- cbind(
              scan_centers,
              scan_weights,
              rep(group_val, scan_n),
              rep(id_val, scan_n),
              rep(group_scans[[group_idx]][scan_idx], scan_n),
              rep(scan_time, scan_n)
            ) |>
              as_tibble(.name_repair = c("minimal"))
          } else if (component_type == "subsample") {
            km <- scan_red$km

            cluster_points <- foreach(cluster_val = seq_len(scan_n)) %do%
              {
                n_obs <- km$size[cluster_val]
                cluster_obs <- scan_data[
                  km$cluster == cluster_val,
                  ,
                  drop = FALSE
                ]

                n_points <- min(n_obs, subsample_size)

                if (n_obs >= subsample_size) {
                  sampled_indices <- sample(n_obs, subsample_size)
                  cluster_points <- cluster_obs[sampled_indices, ]
                  point_weights <- rep(
                    scan_weights[cluster_val] / subsample_size,
                    subsample_size
                  )
                } else {
                  cluster_points <- cluster_obs
                  point_weights <- rep(
                    scan_weights[cluster_val] / n_obs,
                    n_obs
                  )
                }

                cbind(
                  cluster_points,
                  point_weights,
                  rep(group_val, n_points),
                  rep(id_val, n_points),
                  rep(group_scans[[group_idx]][scan_idx], n_points),
                  rep(scan_time, n_points)
                )
              }

            scan_red <- do.call(rbind, cluster_points) |>
              as_tibble(.name_repair = c("minimal"))
          }

          names(scan_red) <- red_names

          p(sprintf(
            "Scan #%i of %i",
            scan_idx,
            length(group_scans[[group_idx]])
          ))

          scan_red
        }
    })

    groups_reduced[[group_idx]] <- do.call(bind_rows, group_reduced_dfs)

    if (case == 1) {
      groups_reduced[[group_idx]] <- groups_reduced[[group_idx]] |>
        mutate(
          X1 = as.numeric(X1),
          X2 = as.numeric(X2),
          weight = as.numeric(weight),
          time_from_bl = as.numeric(time_from_bl),
          partition = 1
        )
    } else if (case == 10) {
      groups_reduced[[group_idx]] <- groups_reduced[[group_idx]] |>
        mutate(
          X1 = as.numeric(X1),
          X2 = as.numeric(X2),
          X3 = as.numeric(X3),
          weight = as.numeric(weight),
          time_from_bl = as.numeric(time_from_bl),
          partition = 1
        )
    }
  }

  list(
    sim_reduced = sim_reduced,
    population_reduced = population_reduced,
    groups_reduced = groups_reduced
  )
}
