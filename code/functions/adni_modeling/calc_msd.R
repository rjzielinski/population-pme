calc_msd <- function(data, projections) {
  require(furrr, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  future_map(
    seq_len(nrow(data)),
    ~ {
      dist_euclidean(
        as.matrix(select(data, x, y, z))[.x, ],
        projections[.x, -1]
      )^2
    },
    .options = furrr_options(seed = TRUE)
  ) %>%
    reduce(c) %>%
    mean()
}
