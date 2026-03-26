silhouette_plot <- function(embeddings, opacity = 0.75, line_width = 2) {
  require(reticulate, quietly = TRUE)

  use_condaenv("lpme")
  pv <- import("pyvista")
  source(
    "~/Documents/brown/research/lpme-project/code/functions/mesh_projection.R"
  )

  time_points <- sort(unique(embeddings$population_embeddings[, 1]))

  palette <- colorRampPalette(hcl.colors(8, palette = "viridis"))
  palette_colors <- palette(151)
  time_colors <- palette_colors[1:101][ceiling((time_points + 1e-10) * 10)]

  lhipp_plot <- pv$Plotter(shape = tuple(2L, 2L))

  for (time_idx in seq_along(time_points)) {
    if (time_points[time_idx] %% 1 == 0) {
      print(time_points[time_idx])

      temp_pop <- population_embeddings[
        population_embeddings[, 1] == time_points[time_idx],
        -1
      ]
      temp_groups <- map(
        group_embeddings,
        ~ .x[.x[, 1] == time_points[time_idx], -1]
      )

      pop_projection3 <- mesh_projection(temp_pop, axis = 3)
      group_projections <- map(temp_groups, ~ mesh_projection(.x, axis = 3))

      lhipp_plot$subplot(0L, 0L)
      lhipp_plot$add_mesh(
        pop_projection3$boundary,
        show_edges = TRUE,
        color = time_colors[time_idx],
        line_width = line_width,
        opacity = opacity,
        label = paste0("Time = ", round(time_points[time_idx], 2))
      )

      if (time_idx == 1) {
        lhipp_plot$add_text(
          "Population",
          position = "upper_left",
          font_size = 12
        )
      }

      if (time_idx == max(which(time_points %% 1 == 0))) {
        lhipp_plot$add_legend(
          bcolor = "white",
          border = TRUE,
          loc = "upper right",
          face = "line"
        )
      }

      lhipp_plot$subplot(0L, 1L)
      lhipp_plot$add_mesh(
        group_projections[[1]]$boundary,
        show_edges = TRUE,
        color = time_colors[time_idx],
        line_width = line_width,
        opacity = opacity,
        label = paste0("Time = ", round(time_points[time_idx], 2))
      )

      if (time_idx == 1) {
        lhipp_plot$add_text(
          group_vals[1],
          position = "upper_left",
          font_size = 12
        )
      }

      if (time_idx == max(which(time_points %% 1 == 0))) {
        lhipp_plot$add_legend(
          bcolor = "white",
          border = TRUE,
          loc = "upper right",
          face = "line"
        )
      }

      lhipp_plot$subplot(1L, 0L)
      lhipp_plot$add_mesh(
        group_projections[[2]]$boundary,
        show_edges = TRUE,
        color = time_colors[time_idx],
        line_width = line_width,
        opacity = opacity,
        label = paste0("Time = ", round(time_points[time_idx], 2))
      )

      if (time_idx == 1) {
        lhipp_plot$add_text(
          group_vals[2],
          position = "upper_left",
          font_size = 12
        )
      }

      if (time_idx == max(which(time_points %% 1 == 0))) {
        lhipp_plot$add_legend(
          bcolor = "white",
          border = TRUE,
          loc = "upper right",
          face = "line"
        )
      }

      lhipp_plot$subplot(1L, 1L)
      lhipp_plot$add_mesh(
        group_projections[[3]]$boundary,
        show_edges = TRUE,
        color = time_colors[time_idx],
        line_width = line_width,
        opacity = opacity,
        label = paste0("Time = ", round(time_points[time_idx], 2))
      )

      if (time_idx == 1) {
        lhipp_plot$add_text(
          group_vals[3],
          position = "upper_left",
          font_size = 12
        )
      }

      if (time_idx == max(which(time_points %% 1 == 0))) {
        lhipp_plot$add_legend(
          bcolor = "white",
          border = TRUE,
          loc = "upper right",
          face = "line"
        )
      }

      lhipp_plot$link_views()
      lhipp_plot$camera_position <- "xy"
    } else {
      next
    }
  }

  lhipp_plot$show()

  lhipp_plot$screenshot(
    "output/additive_mod_silhouette_plot.png"
  )
}
