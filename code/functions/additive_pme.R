additive_pme <- function(
  reduced_data,
  centers,
  init_params,
  weights,
  groups,
  ids,
  scans,
  times,
  partitions,
  lambda = exp(-20:5),
  gamma = exp(-20:5),
  epsilon = 0.01,
  max_iter = 250,
  ssd_ratio_threshold = 5,
  cores = 1,
  verbose = TRUE,
  plot_progress = FALSE
) {
  require(here, quietly = TRUE)
  require(pme, quietly = TRUE)

  source(here("code/functions/fit_adni_additive_model.R"))
  source(here("code/functions/plot_additive_model.R"))
  source(here("code/functions/print_SSD.R"))

  map(
    list.files(here("code/functions/adni_modeling"), full.names = TRUE),
    source
  )

  group_values <- sort(unique(groups))
  id_values <- sort(unique(ids))
  partition_values <- unique(partitions)

  # Fit initial additive model estimate
  additive_model <- list()
  for (partition_idx in seq_along(partition_values)) {
    partition_groups <- groups[partitions == partition_values[partition_idx]]
    partition_ids <- ids[partitions == partition_values[partition_idx]]
    partition_scans <- scans[partitions == partition_values[partition_idx]]
    partition_times <- times[partitions == partition_values[partition_idx]]

    additive_model[[partition_idx]] <- fit_adni_additive_model(
      centers = centers[[partition_idx]],
      params = init_params[[partition_idx]],
      weights = weights[[partition_idx]],
      lambda = lambda,
      gamma = gamma,
      groups = partition_groups,
      ids = partition_ids,
      scans = partition_scans,
      times = partition_times,
      epsilon = epsilon,
      max_iter = max_iter,
      cores = cores
    )
  }

  if (plot_progress == TRUE) {
    plot_additive_model(
      additive_model,
      init_params,
      group_values,
      margin = 0.02
    )
  }

  print("Updating parameters")

  param_list <- update_params(
    additive_model,
    init_params,
    reduced_data,
    centers,
    partition_values,
    group_values,
    id_values
  )

  params <- param_list$params
  center_projections <- param_list$center_projections

  ssd <- calc_ssd(centers, center_projections, partition_values)
  ssd_ratio <- 10 * epsilon

  n <- 1

  ssd_vec <- vector()
  ssd_vec[n] <- ssd

  while (
    (ssd_ratio > epsilon) &
      (ssd_ratio <= ssd_ratio_threshold) &
      (n <= (max_iter - 1))
  ) {
    ssd_prev <- ssd
    additive_model_old <- additive_model
    params_prev <- params

    print("Fitting additive models")

    for (partition_idx in seq_along(partition_values)) {
      partition_groups <- groups[partitions == partition_values[partition_idx]]
      partition_ids <- ids[partitions == partition_values[partition_idx]]
      partition_scans <- scans[partitions == partition_values[partition_idx]]
      partition_times <- times[partitions == partition_values[partition_idx]]

      additive_model[[partition_idx]] <- fit_adni_additive_model(
        centers = centers[[partition_idx]],
        params = params[[partition_idx]],
        weights = weights[[partition_idx]],
        lambda = lambda,
        gamma = gamma,
        groups = partition_groups,
        ids = partition_ids,
        scans = partition_scans,
        times = partition_times,
        epsilon = epsilon,
        max_iter = max_iter,
        cores = cores
      )
    }

    if (plot_progress == TRUE) {
      fig <- plot_additive_model(
        additive_model,
        params,
        group_values,
        margin = 0.02
      )
      print(fig)
    }

    print("Updating parameters")

    param_list <- update_params(
      additive_model,
      params,
      reduced_data,
      centers,
      partition_values,
      group_values,
      id_values
    )

    params <- param_list$params
    center_projections <- param_list$center_projections

    ssd <- calc_ssd(centers, center_projections, partition_values)
    ssd_ratio <- abs(ssd - ssd_prev) / ssd_prev

    n <- n + 1
    ssd_vec[n] <- ssd

    if (verbose == TRUE) {
      print_ssd(ssd, ssd_ratio, n)
    }
  }

  if (plot_progress == TRUE) {
    plot_additive_model(
      additive_model,
      params,
      group_values,
      margin = 0.02
    )
  }

  out_list <- list(
    additive_model = additive_model,
    params = params,
    center_projections = center_projections
  )
}
