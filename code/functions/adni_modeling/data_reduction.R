data_reduction <- function(
  lhipp_surface,
  lhipp_groups,
  lhipp_ids,
  lhipp_scans
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
        scan_data[[partition_idx]] <- scan_df |>
          filter(partition == partition_values[partition_idx]) |>
          select(x, y, z) |>
          as.matrix()

        scan_red[[partition_idx]] <- hdmde(
          scan_data[[partition_idx]],
          50,
          0.05,
          200
        )

        scan_centers[[partition_idx]] <- scan_red[[partition_idx]]$mu
        scan_weights[[partition_idx]] <- scan_red[[partition_idx]]$theta_hat
        scan_n[[partition_idx]] <- nrow(scan_centers[[partition_idx]])

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
        names(scan_red[[partition_idx]]) <- lhipp_red_names
        scan_red[[partition_idx]] <- scan_red[[partition_idx]] |>
          mutate_at(c("x", "y", "z", "weight", "time_from_bl"), as.numeric)
      }

      scan_red_df <- bind_rows(scan_red)
      p(sprintf("Scan #%i of %i", scan_idx, length(lhipp_scans)))
      scan_red_df
    }

  surface_red <- reduce(reduced_dfs, bind_rows)
}
