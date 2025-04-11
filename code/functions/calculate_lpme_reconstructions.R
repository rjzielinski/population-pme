calculate_lpme_reconstructions <- function(lpme_result, data) {
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
  centers <- map(
    1:nrow(params),
    ~ lpme_result$embedding_map(params[.x, ])
  ) |>
    reduce(rbind)

  # within the clusters calculated for each time point, find the nearest
  # cluster center for each observation.
  # adjustments are needed to compare to indices among full set of clusters.
  nearest_center <- map(
    1:nrow(data),
    ~ which.min(
      apply(
        centers[centers[, 1] == data[.x, 1], ],
        1,
        dist_euclidean,
        y = data[.x, ]
      )
    )
  ) |>
    reduce(c)
  center_idx_adjustment <- map(
    1:length(nearest_center),
    ~ nrow(
      centers[centers[, 1] < data[.x, 1], ]
    )
  ) |>
    reduce(c)
  nearest_center <- nearest_center + center_idx_adjustment

  # obtain the parameters associated with each nearest_cluster
  nearest_params <- map(
    1:length(nearest_center),
    ~ params[nearest_center[.x], ]
  ) |>
    reduce(rbind)

  # find d-dimensional projections to the estimated manifold for each observation
  lpme_projections <- map(
    1:nrow(data),
    ~ projection_lpme(
      data[.x, ],
      function(x) lpme_result$embedding_map(x),
      nearest_params[.x, ]
    )
  ) |>
    reduce(rbind)

  # embed projections in D-dimensional space to calculate reconstructions
  lpme_reconstructions <- map(
    1:nrow(data),
    ~ lpme_result$embedding_map(lpme_projections[.x, ])
  ) |>
    reduce(rbind)

  return(
    list(
    reconstructions = lpme_reconstructions,
    projections = lpme_projections
  )
  )
}
