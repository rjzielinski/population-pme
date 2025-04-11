calculate_pme_reconstructions <- function(pme_result, data) {
  #' Calculate PME Reconstructions
  #'
  #' Uses a fitted PME model to estimate reconstructed observations in D-dimensional space
  #'
  #' @param lpme_result The fitted PME model
  #' @param data A numeric matrix of d-dimensional parameters, with time values excluded

  require(purrr)

  # obtain pme run with optimal tuning value and
  # parameter values for the associated cluster centers
  pme_opt <- which.min(pme_result$MSD)
  params <- pme_result$parameterization[[pme_opt]]

  centers <- map(
    seq_len(nrow(params)),
    ~ pme_result$embedding_map(params[.x, ])
  ) |>
    reduce(rbind)

  # calculate the cluster centers nearest to each observation
  nearest_center <- map(
    seq_len(nrow(data)),
    ~ which.min(
      apply(
        centers,
        1,
        dist_euclidean,
        y = data[.x, ]
      )
    )
  ) |>
    reduce(c)

  # obtain the parameters associated with each nearest_cluster
  nearest_params <- map(
    seq_along(nearest_center),
    ~ params[nearest_center[.x], ]
  ) |>
    reduce(rbind)

  # find d-dimensional projections to the estimated manifold for each observation
  pme_projections <- map(
    seq_len(nrow(data)),
    ~ projection_pme(
      data[.x, ],
      function(x) pme_result$embedding_map(x),
      nearest_params[.x, ]
    )
  ) |>
    reduce(rbind)

  # embed projections in D-dimensional space to calculate reconstructions
  pme_reconstructions <- map(
    seq_len(nrow(data)),
    ~ pme_result$embedding_map(pme_projections[.x, ])
  ) |>
    reduce(rbind)

  return(
    list(
      reconstructions = pme_reconstructions, 
      projections = pme_projections
    )
  )
}
