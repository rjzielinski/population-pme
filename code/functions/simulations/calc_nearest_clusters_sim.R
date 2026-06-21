calc_nearest_clusters_sim <- function(
  sim_data,
  reduced_data,
  params,
  cores,
  verbose = FALSE
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

  D <- sim_data |>
    select(contains("X"), -contains("true")) |>
    ncol()
  d <- ncol(params[[1]]) - 1

  sim_data <- sim_data |>
    select(contains("X"), -contains("true"), id, scan, time) |>
    data.table()

  sim_data[, row_idx := .I]

  reduced_data <- data.table(reduced_data)

  scans <- unique(sim_data$scan)

  reduced_ids <- reduced_data[, id]
  reduced_scans <- reduced_data[, scan]
  reduced_times <- reduced_data[, time_from_bl]

  nearest_params <- matrix(
    0,
    nrow = nrow(sim_data),
    ncol = ncol(params[[1]])
  )

  if (verbose == TRUE) {
    p <- progressor(length(scans))
  }
  for (scan_idx in seq_len(length(scans))) {
    scan_data <- sim_data[scan == scans[scan_idx]]
    scan_id <- unique(scan_data$id)

    data_vars <- names(scan_data)[grepl("X", names(scan_data))]
    center_vars <- c(
      "time_from_bl",
      names(reduced_data)[grepl("X", names(reduced_data))]
    )

    scan_centers <- reduced_data[
      scan == scans[scan_idx],
      ..center_vars
    ] |>
      as.matrix()

    scan_params <- params[[1]][reduced_scans == scans[scan_idx], ]

    scan_dist <- dist_betw_matrices(
      as.matrix(scan_data[, ..data_vars]),
      as.matrix(scan_centers[, -1]),
      distance = "rmsd",
      cores = cores
    )

    nearest_param_idx <- rowMins(scan_dist, value = FALSE)

    nearest_params[scan_data$row_idx, ] <- scan_params[
      nearest_param_idx,
    ]

    if (verbose == TRUE) {
      p(sprintf("Scan %d of %d", scan_idx, length(scans)))
    }
  }

  nearest_params
}
