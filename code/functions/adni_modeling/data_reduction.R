data_reduction <- function(
  lhipp_surface,
  lhipp_groups,
  lhipp_ids,
  lhipp_scans,
  min_clusters = 50,
  component_type = "centers",
  subsample_size = 5
) {
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  p <- progressor(steps = length(lhipp_scans))

  partition_values <- unique(lhipp_surface$partition)
  partition_values <- partition_values[order(partition_values)]

  lhipp_red_names <- c(
    "x",
    "y",
    "z",
    "weight",
    "group",
    "id",
    "scan",
    "time_from_bl",
    "partition"
  )

  reduced_dfs <- foreach(
    scan_idx = seq_along(lhipp_scans),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      # print(scan_idx)
      scan_indices <- lhipp_scans == lhipp_scans[scan_idx]

      group_val <- lhipp_surface |>
        filter(image_id == lhipp_scans[scan_idx]) |>
        select(Group) |>
        unique() |>
        unlist()

      id_val <- lhipp_surface |>
        filter(image_id == lhipp_scans[scan_idx]) |>
        select(subid) |>
        unique() |>
        unlist()

      scan_data <- list()
      scan_red <- list()
      scan_centers <- list()
      scan_weights <- list()
      scan_n <- list()

      scan_df <- lhipp_surface |>
        filter(image_id == lhipp_scans[scan_idx])

      scan_time <- unique(scan_df$time_from_bl)

      for (partition_idx in seq_along(partition_values)) {
        # print(partition_idx)

        scan_data[[partition_idx]] <- scan_df |>
          filter(partition == partition_values[partition_idx]) |>
          select(x, y, z) |>
          as.matrix()

        partition_min_clusters <- ifelse(
          nrow(scan_data[[partition_idx]]) <= min_clusters,
          nrow(scan_data[[partition_idx]]) / 2,
          min_clusters
        )

        scan_red[[partition_idx]] <- hdmde(
          scan_data[[partition_idx]],
          partition_min_clusters,
          0.05,
          200
        )

        scan_centers[[partition_idx]] <- scan_red[[partition_idx]]$mu
        scan_weights[[partition_idx]] <- scan_red[[partition_idx]]$theta_hat
        scan_n[[partition_idx]] <- nrow(scan_centers[[partition_idx]])

        if (component_type == "centers") {
          scan_red[[partition_idx]] <- cbind(
            scan_centers[[partition_idx]],
            scan_weights[[partition_idx]],
            rep(group_val, scan_n[[partition_idx]]),
            rep(id_val, scan_n[[partition_idx]]),
            rep(lhipp_scans[scan_idx], scan_n[[partition_idx]]),
            rep(scan_time, scan_n[[partition_idx]]),
            rep(partition_values[partition_idx], scan_n[[partition_idx]])
          ) |>
            as_tibble(.name_repair = c("minimal"))
        } else if (component_type == "subsample") {
          partition_km <- scan_red[[partition_idx]]$km

          cluster_points <- foreach(
            cluster_val = seq_len(scan_n[[partition_idx]])
          ) %do%
            {
              n_obs <- partition_km$size[cluster_val]

              cluster_obs <- scan_data[[partition_idx]][
                partition_km$cluster == cluster_val,
                ,
                drop = FALSE
              ]

              n_points <- min(n_obs, subsample_size)

              if (n_obs >= subsample_size) {
                sampled_indices <- sample(n_obs, subsample_size)
                cluster_points <- cluster_obs[sampled_indices, ]
                point_weights <- rep(
                  scan_weights[[partition_idx]][cluster_val] / subsample_size,
                  subsample_size
                )
              } else {
                cluster_points <- cluster_obs
                point_weights <- rep(
                  scan_weights[[partition_idx]][cluster_val] / n_obs,
                  n_obs
                )
              }

              cbind(
                cluster_points,
                point_weights,
                rep(group_val, n_points),
                rep(id_val, n_points),
                rep(lhipp_scans[scan_idx], n_points),
                rep(scan_time, n_points),
                rep(partition_values[partition_idx], n_points)
              )
            }
          scan_red[[partition_idx]] <- do.call(rbind, cluster_points) |>
            as_tibble(.name_repair = c("minimal"))
        }

        names(scan_red[[partition_idx]]) <- lhipp_red_names
        scan_red[[partition_idx]] <- scan_red[[partition_idx]] |>
          mutate_at(c("x", "y", "z", "weight", "time_from_bl"), as.numeric)
      }

      scan_red_df <- bind_rows(scan_red)
      p(sprintf("Scan #%i of %i", scan_idx, length(lhipp_scans)))
      scan_red_df
    }

  surface_red <- do.call(bind_rows, reduced_dfs)
}
