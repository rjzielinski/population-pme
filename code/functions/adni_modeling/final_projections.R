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

  surface_data[, row_num := .I]

  partition_indices <- match(surface_data$partition, partition_values)
  group_indices <- match(surface_data$Group, group_values)
  id_indices <- match(surface_data$subid, id_values)

  param_grids <- map(
    partition_values,
    ~ as.matrix(additive_model[[.x]]$param_grid)
  )

  data.table::setDTthreads(1)

  with_progress({
    p <- progressor(length(id_values[unique(id_indices)]))
    projection_list <- foreach(
      id_idx = seq_along(id_values[unique(id_indices)]),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        group_idx <- unique(group_indices[id_indices == id_idx])

        id_data <- surface_data[
          id_indices == id_idx,
          .(time_from_bl, x, y, z, partition, row_num)
        ]

        id_out <- list()

        for (partition_val in partition_values) {
          embedding_map <- additive_model[[partition_val]]$embeddings[[
            id_idx
          ]]$embedding_map

          spline_coef_map <- additive_model[[partition_val]]$embeddings[[
            id_idx
          ]]$spline_coef_map

          population_embedding <- additive_model[[
            partition_val
          ]]$population_embedding$embedding_map
          group_embedding <- additive_model[[partition_val]]$group_embeddings[[
            group_idx
          ]]$embedding_map

          part_id_data <- id_data[partition == partition_val]
          part_id_row_nums <- part_id_data$row_num

          part_id_data_red <- part_id_data[, .(time_from_bl, x, y, z)] |>
            as.matrix()

          part_id_nearest_params <- nearest_params[part_id_row_nums, ]

          part_id_params <- map(
            seq_len(nrow(part_id_data)),
            ~ projection_additive_pme(
              part_id_data_red[.x, ],
              spline_coef_map,
              part_id_nearest_params[.x, ],
              param_grids[[partition_val]]
            )
          )

          part_id_params <- do.call(rbind, part_id_params)

          part_id_projections <- apply(part_id_params, 1, embedding_map) |>
            t()

          part_id_pop_projections <- apply(
            part_id_params,
            1,
            population_embedding
          ) |>
            t()

          part_id_group_projections <- apply(
            part_id_params,
            1,
            group_embedding
          ) |>
            t()

          part_id_group_projections <- cbind(
            part_id_params[, 1],
            part_id_pop_projections[, -1] + part_id_group_projections[, -1]
          )

          id_out[[partition_val]] <- cbind(
            part_id_row_nums,
            part_id_params,
            part_id_projections,
            part_id_pop_projections,
            part_id_group_projections
          )
        }

        p(sprintf("ID %d of %d", id_idx, length(id_values)))

        id_out <- do.call(rbind, id_out)
      }
  })

  projection_mat <- do.call(rbind, projection_list)
  projection_mat <- projection_mat[order(projection_mat[, 1]), -1]

  param_mat <- projection_mat[, 1:(d + 1)]
  projections <- projection_mat[, (d + 2):((d + 1) + (D + 1))]

  pop_projections <- projection_mat[,
    ((d + 1) + (D + 1) + 1):((d + 1) + 2 * (D + 1))
  ]

  group_projections <- projection_mat[,
    ((d + 1) + 2 * (D + 1) + 1):((d + 1) + 3 * (D + 1))
  ]

  list(
    params = param_mat,
    projections = projections,
    population_projections = pop_projections,
    group_projections = group_projections
  )
}
