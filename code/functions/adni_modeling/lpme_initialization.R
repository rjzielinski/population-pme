lpme_initialization <- function(data, ids, ...) {
  # fit LPME model to first individual's data
  # assuming that data data.frame object has columns:
  # time_from_bl, x, y, z, subid, partition

  require(dplyr, quietly = TRUE)
  require(doFuture, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(pme, quietly = TRUE)

  partition_values <- unique(data$partition)
  partition_values <- partition_values[order(partition_values)]
  init_data <- list()

  lhipp_init_lpme <- foreach(
    partition_idx = seq_along(partition_values),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      partition_value <- partition_values[partition_idx]
      init_data <- data |>
        filter(
          subid == ids[1],
          partition == partition_value
        ) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()

      lpme(
        init_data,
        ...
      )
    }

  lhipp_init_lpme
}
