calc_ssd <- function(centers, center_projections, partition_values) {
  require(furrr, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(purrr, quietly = TRUE)

  ssd <- list()
  for (partition_idx in seq_along(partition_values)) {
    ssd[[partition_idx]] <- future_map(
      seq_len(nrow(centers[[partition_idx]])),
      ~ dist_euclidean(
        centers[[partition_idx]][.x, ],
        center_projections[[partition_idx]][.x, ]
      )^2
    ) %>%
      reduce(c) %>%
      sum()
  }

  ssd_full <- reduce(ssd, sum)
}
