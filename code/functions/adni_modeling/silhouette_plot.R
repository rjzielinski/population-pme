silhouette_plot <- function(
  additive_model,
  surface_data,
  projections,
  axis = 3,
  filename = NULL,
  opacity = 0.75,
  line_width = 2
) {
  require(reticulate, quietly = TRUE)

  use_condaenv("lpme")
  pv <- import("pyvista")
  source(
    "~/Documents/brown/research/lpme-project/code/functions/mesh_projection.R"
  )

  group_values <- additive_model[[1]]$group_values
  id_values <- additive_model[[1]]$id_values

  n_groups <- length(group_values)

  groups <- surface_data$Group
  group_indices <- match(groups, group_values)

  times <- round(surface_data$time_from_bl, 1)
  time_points <- sort(unique(times))
  time_indices <- match(times, time_points)

  max_time <- max(time_points)

  palette <- colorRampPalette(hcl.colors(8, palette = "viridis"))
  # add 50 to palette to avoid end of spectrum, which is too light to see well
  palette_colors <- palette((10 * max_time) + 50)
  time_colors <- palette_colors[ceiling((time_points + 1e-10) * 10)]

  lhipp_plot <- pv$Plotter(shape = tuple(2L, 2L))

  population_projections <- projections$population_projections

  # Only consider first three groups - AD, CN, and MCI
  # Exclude SMC group due to small sample size
  group_projections_full <- foreach(group_idx = seq_len(3)) %do%
    {
      projections$group_projections[group_indices == group_idx, ]
    }

  group_times <- map(group_projections_full, ~ sort(unique(round(.x[, 1], 1))))

  for (time_idx in seq_along(time_points)) {
    if (time_points[time_idx] %% 1 == 0) {
      print(time_points[time_idx])

      temp_pop <- population_projections[time_indices == time_idx, -1]

      temp_groups <- map(
        group_projections_full,
        ~ .x[round(.x[, 1], 1) == time_points[time_idx], -1]
      )

      pop_projection <- mesh_projection(temp_pop, axis = axis)
      group_projections <- list()

      for (group_idx in seq_along(temp_groups)) {
        if (time_points[time_idx] %in% group_times[[group_idx]]) {
          group_projections[[group_idx]] <- mesh_projection(
            temp_groups[[group_idx]],
            axis = axis
          )
        }
      }

      lhipp_plot$subplot(0L, 0L)
      lhipp_plot$add_mesh(
        pop_projection$boundary,
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

      if (time_points[time_idx] %in% group_times[[1]]) {
        lhipp_plot$subplot(0L, 1L)
        lhipp_plot$add_mesh(
          group_projections[[1]]$boundary,
          show_edges = TRUE,
          color = time_colors[time_idx],
          line_width = line_width,
          opacity = opacity,
          label = paste0("Time = ", round(time_points[time_idx], 2))
        )

        if (
          time_points[time_idx] ==
            min(group_times[[1]][group_times[[1]] %% 1 == 0])
        ) {
          lhipp_plot$add_text(
            group_values[1],
            position = "upper_left",
            font_size = 12
          )
        }

        if (
          time_points[time_idx] ==
            max(group_times[[1]][group_times[[1]] %% 1 == 0])
        ) {
          lhipp_plot$add_legend(
            bcolor = "white",
            border = TRUE,
            loc = "upper right",
            face = "line"
          )
        }
      }

      # if (time_idx == max(which(time_points %% 1 == 0))) {

      if (time_points[time_idx] %in% group_times[[2]]) {
        lhipp_plot$subplot(1L, 0L)
        lhipp_plot$add_mesh(
          group_projections[[2]]$boundary,
          show_edges = TRUE,
          color = time_colors[time_idx],
          line_width = line_width,
          opacity = opacity,
          label = paste0("Time = ", round(time_points[time_idx], 2))
        )

        if (
          time_points[time_idx] ==
            min(group_times[[2]][group_times[[2]] %% 1 == 0])
        ) {
          lhipp_plot$add_text(
            group_values[2],
            position = "upper_left",
            font_size = 12
          )
        }
        if (
          time_points[time_idx] ==
            max(group_times[[2]][group_times[[2]] %% 1 == 0])
        ) {
          lhipp_plot$add_legend(
            bcolor = "white",
            border = TRUE,
            loc = "upper right",
            face = "line"
          )
        }
      }

      if (time_points[time_idx] %in% group_times[[3]]) {
        lhipp_plot$subplot(1L, 1L)
        lhipp_plot$add_mesh(
          group_projections[[3]]$boundary,
          show_edges = TRUE,
          color = time_colors[time_idx],
          line_width = line_width,
          opacity = opacity,
          label = paste0("Time = ", round(time_points[time_idx], 2))
        )

        if (
          time_points[time_idx] ==
            min(group_times[[3]][group_times[[3]] %% 1 == 0])
        ) {
          lhipp_plot$add_text(
            group_values[3],
            position = "upper_left",
            font_size = 12
          )
        }

        if (
          time_points[time_idx] ==
            max(group_times[[3]][group_times[[3]] %% 1 == 0])
        ) {
          lhipp_plot$add_legend(
            bcolor = "white",
            border = TRUE,
            loc = "upper right",
            face = "line"
          )
        }
      }

      if (time_idx == max(which(time_points %% 1 == 0))) {
        lhipp_plot$link_views()
        if (axis == 1) {
          lhipp_plot$camera_position <- "yz"
        } else if (axis == 2) {
          lhipp_plot$camera_position <- "xz"
        } else if (axis == 3) {
          lhipp_plot$camera_position <- "xy"
        }
      }
    } else {
      next
    }
  }

  lhipp_plot$show()

  if (!is.null(filename)) {
    lhipp_plot$screenshot(filename)
  }
}
