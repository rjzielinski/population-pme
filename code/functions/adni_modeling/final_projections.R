final_projections <- function(
  additive_model,
  surface_data,
  reduced_data,
  params,
  partition_values,
  group_values,
  id_values
) {
  require(doFuture, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(here, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  source(here("code/functions/adni_modeling/calc_nearest_clusters.R"))

  nearest_clusters <- calc_nearest_clusters(
    surface_data,
    reduced_data,
    params,
    partition_values
  )

  nearest_params <- nearest_clusters$params

  projections <- foreach(
    row_idx = seq_len(nrow(surface_data)),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      row_id <- surface_data$subid[row_idx]
      row_group <- surface_data$group[row_idx]
      row_partition <- surface_data$partition[row_idx]
      row_point <- surface_data |>
        select(time_from_bl, x, y, z) |>
        slice(row_idx) |>
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

      pop_projection <- additive_model[[partition_idx]]$population_embedding(
        param
      )
      group_projection <- c(
        param[1],
        additive_model[[partition_idx]]$population_embedding(param)[-1] +
          additive_model[[partition_idx]]$group_embedding[[group_idx]](param)[
            -1
          ]
      )

      list(
        param = param,
        full_projection = full_projection,
        population_projection = pop_projection,
        group_projection = group_projection
      )
    }

  full_params <- map(projections, ~ .x$param) |>
    reduce(rbind)
  population_projections <- map(
    projections,
    ~ .x$population_projection
  ) |>
    reduce(rbind)
  group_projections <- map(projections, ~ .x$group_projection) |>
    reduce(rbind)
  projections <- map(projections, ~ .x$full_projection) |>
    reduce(rbind)

  list(
    params = full_params,
    projections = projections,
    population_projections = population_projections,
    group_projections = group_projections
  )
}
