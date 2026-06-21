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
  template,
  lambda = exp(-20:5),
  gamma = exp(-20:5),
  epsilon = 0.025,
  max_iter = 250,
  ssd_ratio_threshold = 5,
  cores = 1,
  verbose = TRUE,
  plot_progress = FALSE
) {
  require(here, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(pracma, quietly = TRUE)

  source(here("code/functions/scale_spherical.R"))
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
    if (template == "sphere") {
      params_rescaled <- init_params[[partition_idx]][,
        -1
      ] |>
        cbind(1) |>
        sph2cart() |>
        cart2sph()

      init_params[[partition_idx]][, -1] <- params_rescaled[, -3]
    }

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
      template = template,
      epsilon = epsilon,
      max_iter = max_iter,
      cores = cores,
      verbose = verbose,
      plot_progress = plot_progress
    )
  }

  if (plot_progress == TRUE) {
    fig <- plot_additive_model(
      additive_model,
      init_params,
      group_values,
      margin = 0.02
    )
    print(fig)
  }

  if (verbose == TRUE) {
    print("Updating parameters")
  }

  param_list <- update_params(
    additive_model,
    init_params,
    reduced_data,
    centers,
    partition_values,
    group_values,
    id_values,
    template,
    verbose = verbose
  )

  params <- param_list$params
  center_projections <- param_list$center_projections

  if (plot_progress == TRUE) {
    fig <- plot_additive_model(
      additive_model,
      params,
      group_values,
      margin = 0.02
    )
    print(fig)
  }

  if (template == "sphere") {
    for (partition_idx in seq_along(partition_values)) {
      params_rescaled <- params[[partition_idx]][, -1] |>
        cbind(1) |>
        sph2cart() |>
        cart2sph()

      params[[partition_idx]][, -1] <- params_rescaled[, -3]
    }
  }

  ssd <- calc_ssd(centers, center_projections, partition_values)
  ssd_ratio <- 10 * epsilon

  if (verbose == TRUE) {
    print_ssd(ssd, NA, 0)
  }

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

    if (verbose == TRUE) {
      print("Fitting additive models")
    }

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
        template = template,
        epsilon = epsilon,
        max_iter = max_iter,
        cores = cores,
        verbose = verbose,
        plot_progress = plot_progress
      )
    }

    if (verbose == TRUE) {
      print("Updating parameters")
    }

    param_list <- update_params(
      additive_model,
      params,
      reduced_data,
      centers,
      partition_values,
      group_values,
      id_values,
      template,
      verbose = verbose
    )

    params <- param_list$params
    center_projections <- param_list$center_projections

    if (template == "sphere") {
      for (partition_idx in seq_along(partition_values)) {
        params_rescaled <- params[[partition_idx]][, -1] |>
          cbind(1) |>
          sph2cart() |>
          cart2sph()

        params[[partition_idx]][, -1] <- params_rescaled[, -3]
      }
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
