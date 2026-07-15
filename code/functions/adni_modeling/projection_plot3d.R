projection_plot3d <- function(
  embeddings,
  filename = NULL,
  axis = 3,
  opacity = 0.75,
  linewidth = 2,
  interval = 0.5
) {
  require(reticulate, quietly = TRUE)

  use_condaenv("lpme")
  pv <- import("pyvista")
  source(
    "~/Documents/brown/research/lpme-project/code/functions/mesh_projection.R"
  )

  population_embeddings <- embeddings$population_embeddings |>
    reduce(rbind)
  time_points <- unique(population_embeddings[, 1])

  population_bounds <- colMinsMaxs(population_embeddings[, -1]) |>
    as.vector()

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

    temp_projections <- population_embeddings[
      population_embeddings[, 1] == selected_time_points[time_idx],
      -1
    ]

    temp_projections_cloud <- pv$PolyData(temp_projections)
    # temp_projections_mesh <- temp_projections_cloud$reconstruct_surface()

    projection_plot$add_mesh(
      temp_projections_cloud,
      scalars = temp_projections[, 3],
      cmap = "viridis",
      show_scalar_bar = FALSE,
      point_size = 8,
      render_points_as_spheres = TRUE,
      smooth_shading = TRUE
    )

    # projection_plot$add_mesh(
    #   temp_projections_mesh,
    #   color = TRUE,
    #   show_edges = TRUE
    # )

    projection_plot$show_grid(
      bounds = population_bounds,
      color = "gray",
      padding = 0.05,
      xtitle = "X",
      ytitle = "Y",
      ztitle = "Z",
      font_size = 10L
    )

    projection_plot$add_text(
      paste0("Time = ", round(selected_time_points[time_idx], 2)),
      position = "upper_left",
      font_size = 12
    )

    # projection_plot$add_text(
    #   paste0("Volume = ", round(temp_projections_mesh$volume, 2)),
    #   position = "lower_left",
    #   font_size = 12
    # )

    projection_plot$camera_position <- list(
      c(6.0, 1.25, 2.0),
      c(0.0, 0.0, 0.0),
      c(0.0, 0.0, 1.0)
    )
  }

  projection_plot$link_views()
  projection_plot$show()

  if (!is.null(filename)) {
    projection_plot$screenshot(filename)
  }
}
