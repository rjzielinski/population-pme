calculate_lpme_reconstructions <- function(lpme_result, data) {
  require(furrr, quietly = TRUE)

  #' Calculate LPME Reconstructions
  #'
  #' Uses a fitted LPME model to estimate reconstructed observations in D-dimensional space
  #'
  #' @param lpme_result The fitted LPME model
  #' @param data A numeric matrix including a column of time values and d-dimensional parameters

  # obtain lpme run with optimal tuning value and
  # parameter values for the associated cluster centers
  lpme_opt <- which.min(lpme_result$msd)
  params <- lpme_result$parameterization_list[[lpme_opt]]

  # convert the cluster center parameters to D dimensions
  centers <- future_map(
    seq_len(nrow(params)),
    ~ lpme_result$embedding_map(params[.x, ]),
    .options = furrr_options(seed = TRUE)
  ) |>
    reduce(rbind)

  # sort by time point
  params <- params[order(params[, 1]), ]
  centers <- centers[order(centers[, 1]), ]

  center_times <- unique(centers[, 1])

  # within the clusters calculated for each time point, find the nearest
  # cluster center for each observation.
  nearest_center <- future_map(
    seq_len(nrow(data)),
    ~ which.min(
      apply(
        centers[
          # consider cluster centers for image with closest time from baseline
          centers[, 1] ==
            center_times[which.min(abs(center_times - data[.x, 1]))],
        ],
        1,
        dist_euclidean,
        y = data[.x, ]
      )
    ),
    .options = furrr_options(seed = TRUE)
  ) |>
    reduce(c)

  # by default, nearest_center gives cluster index within each scan
  # adjustments allow indices to correspond to row in center matrix
  center_idx_adjustment <- future_map(
    seq_along(nearest_center),
    ~ nrow(
      centers[
        centers[, 1] < center_times[which.min(abs(center_times - data[.x, 1]))],
      ]
    ),
    .options = furrr_options(seed = TRUE)
  ) |>
    reduce(c)

  nearest_center <- nearest_center + center_idx_adjustment

  # obtain the parameters associated with each nearest_cluster
  nearest_params <- future_map(
    nearest_center,
    ~ params[.x, ],
    .options = furrr_options(seed = TRUE)
  ) |>
    reduce(rbind)

  # find d-dimensional projections to the estimated manifold for each observation
  lpme_projections <- future_map(
    seq_len(nrow(data)),
    ~ projection_lpme(
      data[.x, ],
      function(x) lpme_result$embedding_map(x),
      nearest_params[.x, ]
    ),
    .options = furrr_options(seed = TRUE)
  ) |>
    reduce(rbind)

  # embed projections in D-dimensional space to calculate reconstructions
  lpme_reconstructions <- future_map(
    seq_len(nrow(data)),
    ~ lpme_result$embedding_map(lpme_projections[.x, ]),
    .options = furrr_options(seed = TRUE)
  ) |>
    reduce(rbind)

  list(reconstructions = lpme_reconstructions, projections = lpme_projections)
}
