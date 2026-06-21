calc_msd <- function(data, projections) {
  require(dplyr, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(Rfast, quietly = TRUE)

  data <- data |>
    select(contains("X_true")) |>
    as.matrix()

  projections <- projections[, -1]

  dist_vec <- rowsums((data - projections)^2, parallel = TRUE)

  mean(dist_vec)
}
