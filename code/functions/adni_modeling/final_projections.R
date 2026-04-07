final_projections <- function(
  additive_model,
  surface_data,
  reduced_data,
  params,
  partition_values,
  group_values,
  id_values,
  d,
  D,
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

  source(here("code/functions/adni_modeling/calc_nearest_clusters.R"))

  surface_data <- as.data.table(surface_data)
  reduced_data <- as.data.table(reduced_data)

  n_rows <- nrow(surface_data)

  nearest_params <- calc_nearest_clusters(
    surface_data,
    reduced_data,
    params,
    partition_values,
    cores
  )

  param_mat <- matrix(0, nrow = n_rows, ncol = d + 1)
  projections <- matrix(0, nrow = n_rows, ncol = D + 1)
  pop_projections <- matrix(0, nrow = n_rows, ncol = D + 1)
  group_projections <- matrix(0, nrow = n_rows, ncol = D + 1)

  partition_indices <- match(surface_data$partition, partition_values)
  group_indices <- match(surface_data$Group, group_values)
  id_indices <- match(surface_data$subid, id_values)

  p <- progressor(n_rows)
  for (row_idx in seq_len(n_rows)) {
    row_point <- surface_data[row_idx, .(time_from_bl, x, y, z)] |>
      unlist()

    embedding_map <- additive_model[[partition_indices[row_idx]]]$embeddings[[
      id_indices[row_idx]
    ]]

    population_embedding <- additive_model[[partition_indices[
      row_idx
    ]]]$population_embedding$embedding_map
    group_embedding <- additive_model[[partition_indices[
      row_idx
    ]]]$group_embeddings[[
      group_indices[row_idx]
    ]]$embedding_map

    param <- projection_lpme(
      row_point,
      embedding_map,
      nearest_params[row_idx, ]
    )

    param_mat[row_idx, ] <- param

    projections[row_idx, ] <- embedding_map(param)
    pop_projection <- population_embedding(param)
    group_projection <- group_embedding(param)

    pop_projections[row_idx, ] <- pop_projection
    group_projections[row_idx, ] <- c(
      param[1],
      pop_projection[-1] + group_projection[-1]
    )

    p(sprintf("Row %d of %d", row_idx, n_rows))
  }

  list(
    params = param_mat,
    projections = projections,
    population_projections = pop_projections,
    group_projections = group_projections
  )
}
