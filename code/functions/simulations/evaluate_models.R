evaluate_models <- function(data, models, case, d, D) {
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)

  true_values <- data$df_true |>
    select(time, contains("X")) |>
    as.matrix()

  observed_data <- data$df_observed |>
    select(time, contains("X")) |>
    as.matrix()

  lpme_error <- map(
    seq_len(nrow(true_values)),
    ~ dist_euclidean(
      true_values[.x, ],
      models$lpme$reconstructions[.x, 1:(1 + D)]
    )^2
  ) |>
    unlist() |>
    mean()

  lpme_time <- models$lpme$fit_time$toc - models$lpme$fit_time$tic

  pme_error <- map(
    seq_len(nrow(true_values)),
    ~ dist_euclidean(
      true_values[.x, ], models$pme$reconstructions[.x, 1:(1 + D)]
    )^2
  ) |>
    unlist() |>
    mean()

  pme_time <- map(
    seq_along(models$pme$fit_time),
    ~ models$pme$fit_time[[.x]]$toc - models$pme$fit_time[[.x]]$tic
  ) |>
    reduce(c) |>
    sum()

  if (!is.null(models$pc$reconstructions)) {
    principal_curve_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$pc$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()

    principal_curve_time <- map(
      seq_along(models$pc$fit_time),
      ~ models$pc$fit_time[[.x]]$toc - models$pc$fit_time[[.x]]$tic
    ) |>
      reduce(c) |>
      sum()
  } else {
    principal_curve_error <- NA
    principal_curve_time <- NA
  }

  if (!is.null(models$lpme_part$reconstructions)) {
    lpme_part_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$lpme_part$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()

    lpme_part_time <- map(
      seq_along(models$lpme_part$fit_time),
      ~ models$lpme_part$fit_time[[.x]]$toc - models$lpme_part$fit_time[[.x]]$tic
    ) |>
      reduce(c) |>
      sum()
  } else {
    lpme_part_error <- NA
    lpme_part_time <- NA
  }

  if (!is.null(models$pme_part$reconstructions)) {
    pme_part_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$pme_part$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()

    pme_part_time1 <- map(
      seq_along(models$pme_part$fit_time[[1]]),
      ~ models$pme_part$fit_time[[1]][[.x]]$toc - models$pme_part$fit_time[[1]][[.x]]$tic
    ) |>
      reduce(c) |>
      sum()
    pme_part_time2 <- map(
      seq_along(models$pme_part$fit_time[[2]]),
      ~ models$pme_part$fit_time[[2]][[.x]]$toc - models$pme_part$fit_time[[2]][[.x]]$tic
    ) |>
      reduce(c) |>
      sum()
    pme_part_time <- pme_part_time1 + pme_part_time2

  } else {
    pme_part_error <- NA
    pme_part_time <- NA
  }

  if (!is.null(models$pc_part$reconstructions)) {
    principal_curve_part_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$pc_part$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()

    principal_curve_part_time1 <- map(
      seq_along(models$pc_part$fit_time[[1]]),
      ~ models$pc_part$fit_time[[1]][[.x]]$toc - models$pc_part$fit_time[[1]][[.x]]$tic
    ) |>
      reduce(c) |>
      sum()
    principal_curve_part_time2 <- map(
      seq_along(models$pc_part$fit_time[[2]]),
      ~ models$pc_part$fit_time[[2]][[.x]]$toc - models$pc_part$fit_time[[2]][[.x]]$tic
    ) |>
      reduce(c) |>
      sum()
    principal_curve_part_time <- principal_curve_part_time1 + principal_curve_part_time2

  } else {
    principal_curve_part_error <- NA
    principal_curve_part_time <- NA
  }

  data_error <- map(
    seq_len(nrow(true_values)),
    ~ dist_euclidean(true_values[.x, ], observed_data[.x, ])^2
  ) |>
    unlist() |>
    mean()

  error_list <- list(
    lpme_error = lpme_error,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    lpme_part_error = lpme_part_error,
    pme_part_error = pme_part_error,
    principal_curve_part_error = principal_curve_part_error,
    data_error = data_error,
    lpme_time = lpme_time,
    pme_time = pme_time,
    principal_curve_time = principal_curve_time,
    lpme_part_time = lpme_part_time,
    pme_part_time = pme_part_time,
    principal_curve_part_time = principal_curve_part_time
  )

  return(error_list)
}
