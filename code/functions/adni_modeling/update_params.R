update_params <- function(
  additive_model,
  prev_params,
  reduced_data,
  centers,
  partition_values,
  id_vals
) {
  require(doFuture, quietly = TRUE)
  require(foreach, quietly = TRUE)
  require(furrr, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE)

  params <- list()
  center_preds <- list()

  for (partition_idx in seq_along(partition_values)) {
    ids <- reduced_data |>
      filter(partition == partition_values[partition_idx]) |>
      pull(id)
    id_indices <- match(ids, id_vals)

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
      embedding_map <- partition_embeddings[[id_indices[row_idx]]]

      param_est <- projection_lpme(
        centers[[partition_idx]][row_idx, ],
        embedding_map,
        prev_params[[partition_idx]][row_idx, ]
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
