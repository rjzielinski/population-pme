final_projections <- function(
  additive_model,
  surface_data,
  reduced_data,
  params,
  partition_values,
  group_values,
  id_values
) {
  require(data.table, quietly = TRUE)
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(here, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  source(here("code/functions/adni_modeling/calc_nearest_clusters.R"))

  nearest_params <- calc_nearest_clusters(
    surface_data,
    reduced_data,
    params,
    partition_values
  )

  surface_data <- as.data.table(surface_data)
  reduced_data <- as.data.table(reduced_data)

  with_progress({
    p <- progressor(nrow(surface_data))
    projections <- foreach(
      row_idx = seq_len(nrow(surface_data)),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        row_id <- surface_data$subid[row_idx]
        row_group <- surface_data$Group[row_idx]
        row_partition <- surface_data$partition[row_idx]

        row_point <- surface_data[row_idx, .(time_from_bl, x, y, z)] |>
          unlist()

        partition_idx <- which(partition_values == row_partition)
        group_idx <- which(group_values == row_group)
        id_idx <- which(id_values == row_id)

        param <- projection_lpme(
          row_point,
          additive_model[[partition_idx]]$embeddings[[id_idx]],
          nearest_params[row_idx, ]
        )

        full_projection <- additive_model[[partition_idx]]$embeddings[[id_idx]](
          param
        )

        pop_projection <- additive_model[[
          partition_idx
        ]]$population_embedding$embedding_map(
          param
        )
        group_projection <- c(
          param[1],
          additive_model[[partition_idx]]$population_embedding$embedding_map(
            param
          )[-1] +
            additive_model[[partition_idx]]$group_embeddings[[
              group_idx
            ]]$embedding_map(
              param
            )[
              -1
            ]
        )

        p()

        c(param, full_projection, pop_projection, group_projection)
      } |>
      reduce(rbind)
  })

  d <- ncol(params[[1]][, -1])
  D <- (ncol(projections) - d) / 3

  full_params <- projections[, 1:d]
  population_projections <- projections[, (d + D + 1):(d + (2 * D))]
  group_projections <- projections[, (d + (2 * D) + 1):(d + (3 * D))]
  projections <- projections[, (d + 1):(d + D)]

  list(
    params = full_params,
    projections = projections,
    population_projections = population_projections,
    group_projections = group_projections
  )
}
