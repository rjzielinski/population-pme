display_test_results <- function(test_results, groups, margin = 0.02) {
  require(dplyr, quietly = TRUE)
  require(plotly, quietly = TRUE)

  population_embeddings <- test_results$embeddings$population_embeddings
  group_embeddings <- test_results$embeddings$group_embeddings
  n_partitions <- length(group_embeddings)
  n_groups <- length(group_embeddings[[1]])
  population_embeddings_full <- reduce(population_embeddings, rbind)
  group_embeddings_full <- list()
  # rejected_full <- reduce(test_results$rejected, c)

  f_param_test_rejected_full <- reduce(test_results$f_test$param_rejected, c)
  f_permute_test_rejected_full <- reduce(
    test_results$f_test$permute_rejected,
    c
  )

  D <- ncol(population_embeddings_full) - 1

  col_names <- c("x", "y", "z")[1:D]

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
      f_param_test_rejected_full,
      f_permute_test_rejected_full
    )
  }

  embeddings_full <- reduce(group_embeddings_full, rbind)
  embeddings_full <- as_tibble(embeddings_full, .name_repair = "minimal")

  names(embeddings_full) <- c(
    "Group",
    "Time",
    col_names,
    paste0("param_rejected_", col_names),
    paste0("permute_rejected_", col_names)
  )

  embeddings_full <- embeddings_full |>
    mutate(
      Group = as.factor(Group),
      Time = as.numeric(Time)
    ) |>
    mutate_at(col_names, as.numeric) |>
    mutate_at(vars(contains("rejected")), as.logical)

  param_group_plots <- list()
  permute_group_plots <- list()
  layout_args <- list()

  param_group_plots_dim1 <- list()
  param_group_plots_dim2 <- list()
  param_group_plots_dim3 <- list()

  permute_group_plots_dim1 <- list()
  permute_group_plots_dim2 <- list()
  permute_group_plots_dim3 <- list()

  for (group_idx in seq_len(n_groups)) {
    group_val <- groups[group_idx]

    scene_name <- paste0("scene", group_idx)
    scene_start <- ((group_idx - 1) / n_groups) + margin
    scene_end <- (group_idx / n_groups) - margin

    layout_args[[scene_name]] <- list(
      domain = list(x = c(scene_start, scene_end))
    )

    if (D == 2) {
      param_group_plots[[group_idx]] <- plot_ly(
        filter(embeddings_full, Group == group_val),
        x = ~Time,
        y = ~x,
        z = ~y,
        color = ~ (param_rejected_x & param_rejected_y),
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

      param_group_plots_dim1[[group_idx]] <- plot_ly(
        filter(embeddings_full, Group == group_val),
        x = ~Time,
        y = ~x,
        z = ~y,
        color = ~param_rejected_x,
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

      param_group_plots_dim2[[group_idx]] <- plot_ly(
        filter(embeddings_full, Group == group_val),
        x = ~Time,
        y = ~x,
        z = ~y,
        color = ~param_rejected_y,
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

      permute_group_plots[[group_idx]] <- plot_ly(
        filter(embeddings_full, Group == group_val),
        x = ~Time,
        y = ~x,
        z = ~y,
        color = ~ (permute_rejected_x & permute_rejected_y),
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

      permute_group_plots_dim1[[group_idx]] <- plot_ly(
        filter(embeddings_full, Group == group_val),
        x = ~Time,
        y = ~x,
        z = ~y,
        color = ~permute_rejected_x,
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

      permute_group_plots_dim2[[group_idx]] <- plot_ly(
        filter(embeddings_full, Group == group_val),
        x = ~Time,
        y = ~x,
        z = ~y,
        color = ~permute_rejected_y,
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

      full_fig <- plot_ly(
        embeddings_full,
        x = ~Time,
        y = ~x,
        z = ~y,
        color = ~Group,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)
      )
    } else if (D == 3) {
      param_group_plots[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~ (param_rejected_x & param_rejected_y & param_rejected_z),
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

      param_group_plots_dim1[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~param_rejected_x,
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

      param_group_plots_dim2[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~param_rejected_y,
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

      param_group_plots_dim3[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~param_rejected_z,
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

      permute_group_plots[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~ (permute_rejected_x &
          permute_rejected_y &
          permute_rejected_z),
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

      permute_group_plots_dim1[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~permute_rejected_x,
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

      permute_group_plots_dim2[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~permute_rejected_y,
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

      permute_group_plots_dim3[[group_idx]] <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~permute_rejected_z,
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

      full_fig <- plot_ly(
        embeddings_full,
        x = ~x,
        y = ~y,
        z = ~z,
        frame = ~Time,
        color = ~Group,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)
      )
    }
  }

  param_group_fig <- subplot(
    param_group_plots,
    nrows = 1,
    margin = margin,
    shareX = TRUE
  ) |>
    layout(layout_args)

  param_group_fig_dim1 <- subplot(
    param_group_plots_dim1,
    nrows = 1,
    margin = margin,
    shareX = TRUE
  ) |>
    layout(layout_args)

  param_group_fig_dim2 <- subplot(
    param_group_plots_dim2,
    nrows = 1,
    margin = margin,
    shareX = TRUE
  ) |>
    layout(layout_args)

  permute_group_fig <- subplot(
    permute_group_plots,
    nrows = 1,
    margin = margin,
    shareX = TRUE
  ) |>
    layout(layout_args)

  permute_group_fig_dim1 <- subplot(
    permute_group_plots_dim1,
    nrows = 1,
    margin = margin,
    shareX = TRUE
  ) |>
    layout(layout_args)

  permute_group_fig_dim2 <- subplot(
    permute_group_plots_dim2,
    nrows = 1,
    margin = margin,
    shareX = TRUE
  ) |>
    layout(layout_args)

  if (D == 2) {
    plots <- list(
      full_plot = full_fig,
      param_group_plot = param_group_fig,
      param_group_plot_dim1 = param_group_fig_dim1,
      param_group_plot_dim2 = param_group_fig_dim2,
      permute_group_plot = permute_group_fig,
      permute_group_plot_dim1 = permute_group_fig_dim1,
      permute_group_plot_dim2 = permute_group_fig_dim2
    )
  } else if (D == 3) {
    param_group_fig_dim3 <- subplot(
      param_group_plots_dim3,
      nrows = 1,
      margin = margin,
      shareX = TRUE
    ) |>
      layout(layout_args)
    permute_group_fig_dim3 <- subplot(
      permute_group_plots_dim3,
      nrows = 1,
      margin = margin,
      shareX = TRUE
    ) |>
      layout(layout_args)

    plots <- list(
      full_plot = full_fig,
      param_group_plot = param_group_fig,
      param_group_plot_dim1 = param_group_fig_dim1,
      param_group_plot_dim2 = param_group_fig_dim2,
      param_group_plot_dim3 = param_group_fig_dim3,
      permute_group_plot = permute_group_fig,
      permute_group_plot_dim1 = permute_group_fig_dim1,
      permute_group_plot_dim2 = permute_group_fig_dim2,
      permute_group_plot_dim3 = permute_group_fig_dim3
    )
  }

  plots
}
