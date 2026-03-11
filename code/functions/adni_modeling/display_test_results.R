display_test_results <- function(test_results, groups, margin = 0.02) {
  require(dplyr, quietly = TRUE)
  require(plotly, quietly = TRUE)

  population_embeddings <- test_results$embeddings$population_embeddings
  group_embeddings <- test_results$embeddings$group_embeddings
  n_partitions <- length(group_embeddings)
  n_groups <- length(group_embeddings[[1]])
  population_embeddings_full <- reduce(population_embeddings, rbind)
  group_embeddings_full <- list()
  rejected_full <- reduce(test_results$rejected, c)

  for (group_idx in seq_len(n_groups)) {
    group_embeddings_full[[group_idx]] <- map(
      group_embeddings,
      ~ .x[[group_idx]]
    ) |>
      reduce(rbind)

    group_embeddings_full[[group_idx]] <- cbind(
      group_embeddings_full[[group_idx]][, 1],
      population_embeddings_full[, -1] +
        group_embeddings_full[[group_idx]][, -1]
    )

    group_embeddings_full[[group_idx]] <- cbind(
      groups[group_idx],
      group_embeddings_full[[group_idx]],
      rejected_full
    )
  }

  embeddings_full <- reduce(group_embeddings_full, rbind)
  embeddings_full <- as_tibble(embeddings_full, .name_repair = "minimal")

  names(embeddings_full) <- c("Group", "Time", "x", "y", "z", "Rejected")

  embeddings_full <- embeddings_full |>
    mutate(
      Group = as.factor(Group),
      Time = as.numeric(Time),
      x = as.numeric(x),
      y = as.numeric(y),
      z = as.numeric(z),
      Rejected = as.logical(Rejected)
    )

  group_plots <- list()
  layout_args <- list()

  for (group_idx in seq_len(n_groups)) {
    group_val <- groups[group_idx]

    scene_name <- paste0("scene", group_idx)
    scene_start <- ((group_idx - 1) / n_groups) + margin
    scene_end <- (group_idx / n_groups) - margin

    layout_args[[scene_name]] <- list(
      domain = list(x = c(scene_start, scene_end))
    )

    group_plots[[group_idx]] <- plot_ly(
      embeddings_full,
      x = ~x,
      y = ~y,
      z = ~z,
      color = ~Rejected,
      frame = ~Time,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3),
      colors = c("#000000", "#FF6666"),
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
