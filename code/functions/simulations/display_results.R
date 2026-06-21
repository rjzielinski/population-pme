display_results <- function(data, models, case, d, D) {
  require(plotly, quietly = TRUE, warn.conflicts = FALSE)
  # TODO: Extend to more simulation cases
  # TODO: Make sure data and true manifold appear
  mat_observed <- as.matrix(data$df_observed)
  mat_true <- as.matrix(data$df_true)
  if (!is.null(models$pc$reconstructions)) {
    if (D == 2) {
      p <- plot_ly(
        x = models$lpme$reconstructions[, 2],
        y = models$lpme$reconstructions[, 1],
        z = models$lpme$reconstructions[, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3, color = "blue"),
        name = "LPME"
      ) |>
        add_markers(
          x = models$pme$reconstructions[, 2],
          y = models$pme$reconstructions[, 1],
          z = models$pme$reconstructions[, 3],
          name = "PME",
          marker = list(size = 3, color = "green")
        ) |>
        add_markers(
          x = models$pc$reconstructions[, 2],
          y = models$pc$reconstructions[, 1],
          z = models$pc$reconstructions[, 3],
          name = "Principal Curve",
          marker = list(size = 3, color = "purple")
        ) |>
        add_markers(
          x = mat_true[, 2],
          y = mat_true[, 1],
          z = mat_true[, 3],
          name = "True Manifold",
          marker = list(size = 3, color = "red")
        ) |>
        add_markers(
          x = mat_observed[, 2],
          y = mat_observed[, 1],
          z = mat_observed[, 3],
          name = "Data",
          opacity = 0.75,
          marker = list(size = 3, color = "black")
        )
    }

    if (D == 3) {
      p <- plot_ly(
        x = models$lpme$reconstructions[, 2],
        y = models$lpme$reconstructions[, 3],
        z = models$lpme$reconstructions[, 4],
        frame = models$lpme$reconstructions[, 1],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3, color = "blue"),
        name = "LPME"
      ) |>
        add_markers(
          x = models$pme$reconstructions[, 2],
          y = models$pme$reconstructions[, 3],
          z = models$pme$reconstructions[, 4],
          frame = models$pme$reconstructions[, 1],
          name = "PME",
          marker = list(size = 3, color = "green")
        ) |>
        add_markers(
          x = models$pc$reconstructions[, 2],
          y = models$pc$reconstructions[, 3],
          z = models$pc$reconstructions[, 4],
          frame = models$pc$reconstructions[, 1],
          name = "Principal Curve",
          marker = list(size = 3, color = "purple")
        ) |>
        add_markers(
          x = mat_true[, 2],
          y = mat_true[, 3],
          z = mat_true[, 4],
          frame = mat_true[, 1],
          name = "True Manifold",
          marker = list(size = 3, color = "red")
        ) |>
        add_markers(
          x = mat_observed[, 2],
          y = mat_observed[, 3],
          z = mat_observed[, 4],
          frame = mat_observed[, 1],
          name = "Data",
          opacity = 0.3,
          marker = list(size = 3, color = "black")
        )
    }
  } else {
    if (D == 2) {
      p <- plot_ly(
        x = models$lpme$reconstructions[, 2],
        y = models$lpme$reconstructions[, 1],
        z = models$lpme$reconstructions[, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3, color = "blue"),
        name = "LPME"
      ) |>
        add_markers(
          x = models$pme$reconstructions[, 2],
          y = models$pme$reconstructions[, 1],
          z = models$pme$reconstructions[, 3],
          name = "PME",
          marker = list(size = 3, color = "green")
        ) |>
        add_markers(
          x = mat_true[, 2],
          y = mat_true[, 1],
          z = mat_true[, 3],
          name = "True Manifold",
          marker = list(size = 3, color = "red")
        ) |>
        add_markers(
          x = mat_observed[, 2],
          y = mat_observed[, 1],
          z = mat_observed[, 3],
          name = "Data",
          opacity = 0.2,
          marker = list(size = 3, color = "black")
        )
    }

    if (D == 3) {
      p <- plot_ly(
        x = models$lpme$reconstructions[, 2],
        y = models$lpme$reconstructions[, 3],
        z = models$lpme$reconstructions[, 4],
        frame = models$lpme$reconstructions[, 1],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3, color = "blue"),
        name = "LPME"
      ) |>
        add_markers(
          x = models$pme$reconstructions[, 2],
          y = models$pme$reconstructions[, 3],
          z = models$pme$reconstructions[, 4],
          frame = models$pme$reconstructions[, 1],
          name = "PME",
          marker = list(size = 3, color = "green")
        ) |>
        add_markers(
          x = mat_true[, 2],
          y = mat_true[, 3],
          z = mat_true[, 4],
          frame = mat_true[, 1],
          name = "True Manifold",
          marker = list(size = 3, color = "red")
        ) |>
        add_markers(
          x = mat_observed[, 2],
          y = mat_observed[, 3],
          z = mat_observed[, 4],
          frame = mat_observed[, 1],
          name = "Data",
          marker = list(size = 3, color = "black")
        )
    }
  }

  return(p)
}
