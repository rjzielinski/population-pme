prepare_data <- function(surface_data, reduced_data) {
  partition_values <- unique(surface_data$partition)
  partition_values <- partition_values[order(partition_values)]

  group_values <- unique(reduced_data$group)
  group_values <- group_values[order(group_values)]

  id_values <- unique(reduced_data$id)
  id_values <- id_values[order(id_values)]

  lhipp_x <- list()
  lhipp_centers <- list()
  lhipp_weights <- list()

  for (partition_idx in seq_along(partition_values)) {
    lhipp_x[[partition_idx]] <- surface_data |>
      filter(partition == partition_values[partition_idx]) |>
      select(time_from_bl, x, y, z) |>
      as.matrix()

    lhipp_centers[[partition_idx]] <- reduced_data |>
      filter(partition == partition_values[partition_idx]) |>
      select(time_from_bl, x, y, z) |>
      as.matrix()

    lhipp_weights[[partition_idx]] <- reduced_data |>
      filter(partition == partition_values[partition_idx]) |>
      select(weight) |>
      reduce(c)
  }

  out_list <- list(
    lhipp_x = lhipp_x,
    lhipp_centers = lhipp_centers,
    lhipp_weights = lhipp_weights,
    partition_values = partition_values,
    group_values = group_values,
    id_values = id_values
  )
  list2env(out_list, envir = .GlobalEnv)
}
