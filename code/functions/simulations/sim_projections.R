sim_projections <- function(
  additive_model,
  sim_data,
  reduced_data,
  params,
  group_values,
  id_values,
  d,
  D,
  template,
  cores,
  verbose = FALSE
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

  group_indices <- match(sim_data$group, group_values)
  id_indices <- match(sim_data$id, id_values)

  param_grid <- as.matrix(additive_model[[1]]$param_grid)

  data.table::setDTthreads(1)

  with_progress({
    if (verbose == TRUE) {
      p <- progressor(length(id_values[unique(id_indices)]))
    }
    projection_list <- foreach(
      id_idx = seq_along(id_values[unique(id_indices)]),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        group_idx <- unique(group_indices[id_indices == id_idx])

        data_vars <- c(
          "time",
          names(sim_data)[grepl("X", names(sim_data))],
          "row_num"
        )
        data_vars <- data_vars[-which(grepl("true", data_vars))]

        reduced_vars <- data_vars[-which(grepl("row", data_vars))]

        id_data <- sim_data[
          id_indices == id_idx,
          ..data_vars
        ]

        id_out <- list()

        embedding_map <- additive_model[[1]]$embeddings[[
          id_idx
        ]]$embedding_map

        spline_coef_map <- additive_model[[1]]$embeddings[[
          id_idx
        ]]$spline_coef_map

        population_embedding <- additive_model[[
          1
        ]]$population_embedding$embedding_map
        group_embedding <- additive_model[[1]]$group_embeddings[[
          group_idx
        ]]$embedding_map

        id_row_nums <- id_data$row_num

        id_data_red <- id_data[, ..reduced_vars] |>
          as.matrix()

        id_nearest_params <- nearest_params[id_row_nums, ]

        id_params <- map(
          seq_len(nrow(id_data)),
          ~ projection_additive_pme(
            id_data_red[.x, ],
            spline_coef_map,
            id_nearest_params[.x, ],
            param_grid,
            template
          )
        )

        id_params <- do.call(rbind, id_params)

        id_projections <- apply(id_params, 1, embedding_map) |>
          t()

        id_pop_projections <- apply(id_params, 1, population_embedding) |>
          t()

        id_group_projections <- apply(id_params, 1, group_embedding) |>
          t()

        id_group_projections <- cbind(
          id_params[, 1],
          id_pop_projections[, -1] + id_group_projections[, -1]
        )

        if (verbose == TRUE) {
          p(sprintf("ID %d of %d", id_idx, length(id_values)))
        }

        id_out <- cbind(
          id_row_nums,
          id_params,
          id_projections,
          id_pop_projections,
          id_group_projections
        )
      }
  })

  projection_mat <- do.call(rbind, projection_list)
  projection_mat <- projection_mat[order(projection_mat[, 1]), -1]

  param_dim <- d + 1
  param_mat <- projection_mat[, 1:param_dim]

  proj_dim_start <- d + 2
  proj_dim_end <- (d + 1) + (D + 1)
  projections <- projection_mat[, proj_dim_start:proj_dim_end]

  pop_proj_dim_start <- (d + 1) + (D + 1) + 1
  pop_proj_dim_end <- (d + 1) + 2 * (D + 1)
  pop_projections <- projection_mat[, pop_proj_dim_start:pop_proj_dim_end]

  group_proj_dim_start <- (d + 1) + 2 * (D + 1) + 1
  group_proj_dim_end <- (d + 1) + 3 * (D + 1)
  group_projections <- projection_mat[, group_proj_dim_start:group_proj_dim_end]

  list(
    params = param_mat,
    projections = projections,
    population_projections = pop_projections,
    group_projections = group_projections
  )
}
