calc_nearest_clusters <- function(
  surface_data,
  reduced_data,
  params,
  partition_values
) {
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  p <- progressor(nrow(surface_data))
  nearest_clusters <- foreach(row_idx = seq_len(nrow(surface_data))) %dofuture%
    {
      row_id <- surface_data$subid[row_idx]
      row_scan <- surface_data$image_id[row_idx]
      row_time <- surface_data$time_from_bl[row_idx]
      row_partition <- surface_data$partition[row_idx]

      id_vals <- reduced_data |>
        filter(partition == row_partition) |>
        pull(id)

      scan_vals <- reduced_data |>
        filter(partition == row_partition) |>
        pull(scan)

      time_vals <- reduced_data |>
        filter(partition == row_partition) |>
        pull(time_from_bl)

      row_centers <- reduced_data |>
        filter(
          id == row_id,
          scan == row_scan,
          time_from_bl == row_time,
          partition == row_partition
        ) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()

      row_params <- params[[which(partition_values == row_partition)]][
        id_vals == row_id &
          scan_vals == row_scan &
          time_vals == row_time,
      ]

      surface_loc <- surface_data |>
        select(x, y, z) |>
        slice(row_idx) |>
        unlist()

      center_dists <- map(
        seq_len(nrow(row_centers)),
        ~ dist_euclidean(
          surface_loc,
          row_centers[.x, -1]
        )
      ) |>
        reduce(c)

      center_idx <- which.min(center_dists)
      center_loc <- row_centers[center_idx, ]
      center_param <- row_params[center_idx, ]

      p()

      list(
        center = center_loc,
        param = center_param
      )
    }

  nearest_params <- future_map(nearest_clusters, ~ .x$param) |>
    reduce(rbind)
  nearest_centers <- future_map(nearest_clusters, ~ .x$center) |>
    reduce(rbind)

  list(
    centers = nearest_centers,
    params = nearest_params
  )
}
