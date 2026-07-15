group_comp_silhouette_plot <- function(
  embeddings,
  filename = NULL,
  axis = 3,
  opacity = 0.75,
  line_width = 2,
  interval = 0.5
) {
  require(reticulate, quietly = TRUE)

  use_condaenv("lpme")
  pv <- import("pyvista")
  source(
    "~/Documents/brown/research/lpme-project/code/functions/mesh_projection.R"
  )

  n_partitions <- length(embeddings$population_embeddings)
  n_groups <- length(embeddings$group_embeddings[[1]])
  group_names <- sort(unique(embeddings$group_assignments))

  population_embeddings <- embeddings$population_embeddings |>
    reduce(rbind)

  group_embeddings <- list()
  for (group_idx in seq_len(n_groups)) {
    temp_group_embeddings <- map(
      seq_len(n_partitions),
      ~ embeddings$group_embeddings[[.x]][[group_idx]]
    ) |>
      reduce(rbind)
    group_embeddings[[group_idx]] <- cbind(
      population_embeddings[, 1],
      population_embeddings[, -1] + temp_group_embeddings[, -1]
    )
  }

  full_embeddings <- append(
    list(population_embeddings),
    group_embeddings
  )

  embedding_colors <- c("black", "purple", "green", "blue")
  embedding_names <- c("Population", group_names)

  time_points <- unique(population_embeddings[, 1])

  target_time_points <- seq(
    from = min(time_points),
    to = max(time_points),
    by = interval
  )

  selected_time_points <- map(
    target_time_points,
    ~ time_points[which.min(abs(time_points - .x))]
  ) |>
    reduce(c)

  n_cols <- ceiling(length(selected_time_points) / 2)
  projection_plot <- pv$Plotter(
    shape = tuple(2L, as.integer(n_cols)),
    window_size = c(1500L, 1000L)
  )

  for (time_idx in seq_along(selected_time_points)) {
    plot_col <- (time_idx - 1) %% n_cols
    if (time_idx / (2 * n_cols) <= 0.5) {
      projection_plot$subplot(0L, as.integer(plot_col))
    } else {
      projection_plot$subplot(1L, as.integer(plot_col))
    }

    for (embedding_idx in seq_along(full_embeddings)) {
      temp_embeddings <- full_embeddings[[embedding_idx]]
      time_embeddings <- temp_embeddings[
        temp_embeddings[, 1] == selected_time_points[time_idx],
        -1
      ]

      time_projections <- mesh_projection(time_embeddings, axis = axis)

      projection_plot$add_mesh(
        time_projections$boundary,
        show_edges = TRUE,
        color = embedding_colors[embedding_idx],
        line_width = line_width,
        opacity = opacity,
        label = embedding_names[embedding_idx]
      )
    }

    projection_plot$show_grid(
      color = "gray",
      xtitle = "X",
      ytitle = "Y",
      ztitle = "Z"
    )

    projection_plot$add_text(
      paste0("Time = ", round(selected_time_points[time_idx], 2)),
      position = "upper_left",
      font_size = 10
    )

    if (length(selected_time_points) %% 2 == 0) {
      projection_plot$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    if (axis == 1) {
      projection_plot$camera_position <- "yz"
    } else if (axis == 2) {
      projection_plot$camera_position <- "xz"
    } else if (axis == 3) {
      projection_plot$camera_position <- "xy"
    }
  }

  embedding_colors <- c("black", "purple", "green", "blue")
  embedding_names <- c("Population", group_names)

  legend_entries <- list(
    list("Population", "black"),
    list("AD", "purple"),
    list("CN", "green"),
    list("MCI", "blue")
  )

  # 3. Add the legend, scaling it up to act as a standalone figure element

  if (length(selected_time_points) %% 2 == 1) {
    projection_plot$subplot(1L, as.integer(n_cols - 1))
    projection_plot$add_legend(
      labels = legend_entries,
      bcolor = "white",
      border = FALSE,
      size = c(0.5, 0.35),
      loc = "center left",
      face = "line"
    )
  }

  projection_plot$link_views()
  projection_plot$camera$zoom(1.2)
  projection_plot$show()

  if (!is.null(filename)) {
    projection_plot$screenshot(filename)
  }
}
