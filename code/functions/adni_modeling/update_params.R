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

    p <- progressor(nrow(centers[[partition_idx]]))
    params[[partition_idx]] <- future_map(
      seq_len(nrow(centers[[partition_idx]])),
      ~ {
        p()
        projection_lpme(
          centers[[partition_idx]][.x, ],
          additive_model[[partition_idx]]$embeddings[[which(
            id_vals == ids[.x]
          )]],
          prev_params[[partition_idx]][.x, ]
        )
      },
      .options = furrr_options(seed = TRUE)
    ) |>
      reduce(rbind)

    p <- progressor(nrow(params[[partition_idx]]))
    center_preds[[partition_idx]] <- future_map(
      seq_len(nrow(params[[partition_idx]])),
      ~ {
        p()
        additive_model[[partition_idx]]$embeddings[[which(
          id_vals == ids[.x]
        )]](params[[partition_idx]][.x, ])
      },
      .options = furrr_options(seed = TRUE)
    ) |>
      reduce(rbind)
  }

  list(
    params = params,
    center_projections = center_preds
  )
}
