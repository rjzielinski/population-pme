update_params <- function(
  additive_model,
  prev_params,
  reduced_data,
  centers,
  partition_values,
  group_values,
  id_values
) {
  require(doFuture, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(here, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  source(here("code/functions/projection_additive_pme.R"))

  d <- ncol(prev_params[[1]]) - 1

  params <- list()
  center_preds <- list()

  for (partition_idx in seq_along(partition_values)) {
    param_grid <- additive_model[[partition_idx]]$param_grid |>
      as.matrix()
    n_knots <- nrow(param_grid)

    partition_data <- reduced_data |>
      filter(partition == partition_values[partition_idx])
    groups <- partition_data$group
    ids <- partition_data$id

    id_indices <- match(ids, id_values)
    group_indices <- match(groups, group_values)

    n_rows <- nrow(prev_params[[partition_idx]])

    param_mat <- matrix(
      0,
      nrow = n_rows,
      ncol = ncol(prev_params[[partition_idx]])
    )
    center_pred_mat <- matrix(
      0,
      nrow = n_rows,
      ncol = ncol(centers[[partition_idx]])
    )

    partition_embeddings <- additive_model[[partition_idx]]$embeddings

    p <- progressor(nrow(centers[[partition_idx]]))
    for (row_idx in seq_len(n_rows)) {
      center_val <- centers[[partition_idx]][row_idx, ]
      time_val <- center_val[1]

      embedding_map <- partition_embeddings[[id_indices[row_idx]]]$embedding_map
      spline_coef_map <- partition_embeddings[[id_indices[
        row_idx
      ]]]$spline_coef_map

      param_est <- projection_additive_pme(
        center_val,
        spline_coef_map,
        prev_params[[partition_idx]][row_idx, ],
        param_grid
      )

      param_mat[row_idx, ] <- param_est
      center_pred_mat[row_idx, ] <- embedding_map(param_mat[row_idx, ])
      p()
    }

    params[[partition_idx]] <- param_mat
    center_preds[[partition_idx]] <- center_pred_mat
  }

  list(
    params = params,
    center_projections = center_preds
  )
}
