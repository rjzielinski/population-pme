calc_nearest_clusters <- function(
  surface_data,
  reduced_data,
  params,
  partition_values
) {
  require(data.table, quietly = TRUE)
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(dtplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  surface_data <- data.table(surface_data)
  reduced_data <- data.table(reduced_data)

  p <- progressor(nrow(surface_data))

  nearest_params <- future_map(
    seq_len(nrow(surface_data)),
    function(row_idx) {
      row_id <- surface_data$subid[row_idx]
      row_scan <- surface_data$image_id[row_idx]
      row_time <- surface_data$time_from_bl[row_idx]
      row_partition <- surface_data$partition[row_idx]

      id_vals <- reduced_data[partition == row_partition, id]
      scan_vals <- reduced_data[partition == row_partition, scan]
      time_vals <- reduced_data[partition == row_partition, time_from_bl]

      row_centers <- reduced_data[
        (id == row_id) &
          (scan == row_scan) &
          (time_from_bl == row_time) &
          (partition == row_partition),
        .(time_from_bl, x, y, z)
      ] |>
        as.matrix()

      row_params <- params[[which(partition_values == row_partition)]][
        id_vals == row_id &
          scan_vals == row_scan &
          time_vals == row_time,
      ]

      surface_loc <- surface_data[row_idx, .(x, y, z)] |>
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
      p(sprintf("Row %d", row_idx))

      row_params[center_idx, ]
    }
  ) |>
    reduce(rbind)

  nearest_params
}
