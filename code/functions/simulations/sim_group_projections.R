sim_group_projections <- function(
  additive_model,
  sim_data,
  reduced_data,
  params,
  group_idx,
  d,
  D,
  template,
  cores
) {
  require(data.table, quietly = TRUE)
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(here, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  source(here("code/functions/simulations/calc_nearest_clusters_sim.R"))

  sim_data <- as.data.table(sim_data)
  reduced_data <- as.data.table(reduced_data)

  n_rows <- nrow(sim_data)

  nearest_params <- calc_nearest_clusters_sim(
    sim_data,
    reduced_data,
    params,
    cores
  )

  sim_data[, row_num := .I]

  param_grid <- as.matrix(additive_model[[1]]$param_grid)

  data.table::setDTthreads(1)

  data_vars <- c(
    "time",
    names(sim_data)[grepl("X", names(sim_data))],
    "row_num"
  )
  data_vars <- data_vars[-which(grepl("true", data_vars))]
  reduced_vars <- data_vars[-which(grepl("row", data_vars))]

  data <- sim_data[, ..data_vars]

  proj_out <- list()
  population_embedding <- additive_model[[1]]$population_embedding$embedding_map
  pop_spline_coef_map <- additive_model[[
    1
  ]]$population_embedding$spline_coef_map

  group_embedding <- additive_model[[1]]$group_embeddings[[
    group_idx
  ]]$embedding_map
  group_spline_coef_map <- additive_model[[1]]$group_embeddings[[
    group_idx
  ]]$spline_coef_map

  row_nums <- data$row_num
  data_red <- data[, ..reduced_vars] |>
    as.matrix()

  additive_mod_params <- map(
    seq_len(nrow(data)),
    ~ projection_additive_pme(
      data_red[.x, ],
      function(x) pop_spline_coef_map(x) + group_spline_coef_map(x),
      nearest_params[.x, ],
      param_grid,
      template
    )
  )
  param_mat <- do.call(rbind, additive_mod_params)

  population_projections <- apply(param_mat, 1, population_embedding) |>
    t()
  group_projections <- apply(param_mat, 1, group_embedding) |>
    t()

  projection_mat <- cbind(
    population_projections[, 1],
    population_projections[, -1] + group_projections[, -1]
  )

  list(
    params = param_mat,
    projections = projection_mat
  )
}
