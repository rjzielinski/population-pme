plot_additive_model <- function(
  additive_model,
  params,
  group_vals,
  margin = 0.02
) {
  require(plotly, quietly = TRUE)
  require(purrr, quietly = TRUE)

  n_partitions <- length(additive_model)
  n_groups <- length(group_vals)

  pop_partition_embedding <- list()
  group_embeddings <- lapply(seq_len(n_groups), function(x) list())

  for (partition_idx in seq_len(n_partitions)) {
    partition_params <- params[[partition_idx]]
    pop_embedding_map <- additive_model[[
      partition_idx
    ]]$population_embedding$embedding_map

    pop_embedding <- map(
      seq_len(nrow(params[[partition_idx]])),
      ~ pop_embedding_map(partition_params[.x, ])
    )
    pop_embedding <- do.call(rbind, pop_embedding)
    pop_partition_embedding[[partition_idx]] <- pop_embedding

    for (group_idx in seq_len(n_groups)) {
      group_embedding_map <- additive_model[[partition_idx]]$group_embeddings[[
        group_idx
      ]]$embedding_map

      group_embedding <- map(
        seq_len(nrow(params[[partition_idx]])),
        ~ group_embedding_map(partition_params[.x, ])
      )
      group_embedding <- do.call(rbind, group_embedding)

      group_embeddings[[group_idx]][[partition_idx]] <- cbind(
        pop_embedding[, 1],
        pop_embedding[, -1] + group_embedding[, -1]
      )
    }
  }

  population_embeddings <- do.call(rbind, pop_partition_embedding)
  group_embeddings <- lapply(group_embeddings, function(partition_list) {
    do.call(rbind, partition_list)
  })

  group_plots <- list()
  layout_args <- list()

  for (group_idx in seq_len(n_groups)) {
    group_val <- group_vals[group_idx]

    scene_name <- paste0("scene", group_idx)
    scene_start <- ((group_idx - 1) / n_groups) + margin
    scene_end <- (group_idx / n_groups) - margin

    layout_args[[scene_name]] <- list(
      domain = list(x = c(scene_start, scene_end))
    )

    group_plots[[group_idx]] <- plot_ly(
      x = group_embeddings[[group_idx]][, 2],
      y = group_embeddings[[group_idx]][, 3],
      z = group_embeddings[[group_idx]][, 4],
      frame = group_embeddings[[group_idx]][, 1],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3, color = "black"),
      name = group_val,
      scene = paste0("scene", group_idx)
    ) |>
      add_annotations(
        x = (scene_start + scene_end) / 2,
        y = 1,
        text = group_val,
        xanchor = "center",
        yanchor = "top",
        showarrow = FALSE,
        font = list(size = 15)
      )
  }

  fig <- subplot(group_plots, nrows = 1, margin = margin, shareX = TRUE)

  fig <- fig |>
    layout(layout_args)

  fig
}
