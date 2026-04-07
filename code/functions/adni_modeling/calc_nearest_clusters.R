calc_nearest_clusters <- function(
  surface_data,
  reduced_data,
  params,
  partition_values,
  cores
) {
  require(data.table, quietly = TRUE)
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(dtplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(lineup2, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(Rfast, quietly = TRUE)

  surface_data <- data.table(surface_data)
  surface_data <- surface_data[, .(
    x,
    y,
    z,
    subid,
    image_id,
    time_from_bl,
    partition
  )]

  surface_data[, row_idx := .I]

  reduced_data <- data.table(reduced_data)

  scans <- unique(surface_data$image_id)

  partition_ids <- list()
  partition_scans <- list()
  partition_times <- list()

  for (partition_idx in seq_along(partition_values)) {
    partition_ids[[partition_idx]] <- reduced_data[
      partition == partition_idx,
      id
    ]
    partition_scans[[partition_idx]] <- reduced_data[
      partition == partition_idx,
      scan
    ]
    partition_times[[partition_idx]] <- reduced_data[
      partition == partition_idx,
      time_from_bl
    ]
  }

  nearest_params <- matrix(
    0,
    nrow = nrow(surface_data),
    ncol = ncol(params[[1]])
  )

  p <- progressor(length(scans))
  for (scan_idx in seq_len(length(scans))) {
    scan_data <- surface_data[image_id == scans[scan_idx]]
    scan_partition_vals <- unique(scan_data$partition)
    scan_id <- unique(scan_data$subid)

    for (partition_val in scan_partition_vals) {
      scan_part_data <- scan_data[partition == partition_val]

      scan_part_centers <- reduced_data[
        (scan == scans[scan_idx]) & (partition == partition_val),
        .(time_from_bl, x, y, z)
      ] |>
        as.matrix()

      scan_part_params <- params[[which(partition_values == partition_val)]][
        (partition_scans[[partition_val]] == scans[scan_idx]),
      ]

      scan_part_dist <- dist_betw_matrices(
        as.matrix(scan_part_data[, .(x, y, z)]),
        as.matrix(scan_part_centers[, -1]),
        distance = "rmsd",
        cores = cores
      )

      nearest_param_idx <- rowMins(scan_part_dist, value = FALSE)

      nearest_params[scan_part_data$row_idx, ] <- scan_part_params[
        nearest_param_idx,
      ]
    }

    p(sprintf("Scan %d of %d", scan_idx, length(scans)))
  }

  nearest_params
}
