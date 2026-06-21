fit_models <- function(data, case, d, D, partition = FALSE) {
  ##### MANIFOLD LEARNING #####

  require(pme, quietly = TRUE, warn.conflicts = FALSE)
  require(princurve, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  require(tictoc, quietly = TRUE, warn.conflicts = FALSE)

  mat <- as.matrix(data$df_observed)

  if (case %in% c(3, 8, 9)) {
    partition <- TRUE
  }

  # fit LPME algorithm and calculate reconstructions
  tic()
  lpme_result <- lpme(mat, d, verbose = FALSE, print_plots = FALSE)
  lpme_time <- toc()
  lpme_reconstructions <- calculate_lpme_reconstructions(
    lpme_result,
    mat
  )

  # fit PME and appropriate principal curve/surface algorithms and
  # calculate reconstructions for each
  pme_result_list <- list()
  pc_result_list <- list()

  pme_reconstruction_list <- list()
  pc_reconstruction_list <- list()

  pme_time_list <- list()
  pc_time_list <- list()

  smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")

  time_values <- unique(data$df$time)
  for (time_idx in seq_along(time_values)) {
    temp_data <- mat[mat[, 1] == time_values[time_idx], -1]
    tic()
    pme_result_list[[time_idx]] <- pme(temp_data, d = d)
    pme_time_list[[time_idx]] <- toc()
    pme_reconstruction_list[[time_idx]] <- calculate_pme_reconstructions(
      pme_result_list[[time_idx]],
      temp_data
    )
    pme_reconstruction_list[[time_idx]] <- cbind(
      time_values[time_idx],
      pme_reconstruction_list[[time_idx]]
    )
    if (d == 1) {
      principal_curves <- list()
      pc_times <- list()
      pc_error <- vector()
      for (smoother_idx in seq_along(smoothing_options)) {
        tic()
        principal_curves[[smoother_idx]] <- principal_curve(
          temp_data,
          smoother = smoothing_options[[smoother_idx]]
        )
        pc_times[[smoother_idx]] <- toc()
        pc_error[smoother_idx] <- principal_curves[[smoother_idx]]$dist
      }
      opt_principal_curve <- which.min(pc_error)
      pc_result_list[[time_idx]] <- principal_curves[[opt_principal_curve]]
      pc_reconstruction_list[[time_idx]] <- cbind(
        time_values[time_idx],
        pc_result_list[[time_idx]]$s
      )
      pc_time_list[[time_idx]] <- pc_times[[opt_principal_curve]]
    } else if (d == 2) {
      if (dim(temp_data)[2] == 3) {
        tic()
        principal_surface <- prinSurf(temp_data)
        pc_time_list[[time_idx]] <- toc()
        surface_mse <- map(
          seq_along(principal_surface),
          ~ principal_surface[[.x]]$MSE
        ) |>
          unlist()
        opt_surface <- which.min(surface_mse)
        # using opt_surface + 2 as index because first two list entries are NULL
        pc_result_list[[time_idx]] <- principal_surface[[opt_surface + 2]]
        pc_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface[[opt_surface + 2]]$PS
        )
      } else if (dim(temp_data)[2] > 3) {
        # As coded, the principal surface function assumes that D = 3
        # This is not the case when considering augmented data
        pc_result_list[[time_idx]] <- NULL
        pc_reconstruction_list[[time_idx]] <- NULL
        pc_time_list[[time_idx]] <- NULL
      }
    }
  }

  pme_reconstructions <- reduce(pme_reconstruction_list, rbind)
  if (length(pc_reconstruction_list) > 0) {
    pc_reconstructions <- reduce(pc_reconstruction_list, rbind)
  } else {
    pc_reconstructions <- NULL
  }

  lpme_out <- list(
    lpme = lpme_result,
    reconstructions = lpme_reconstructions,
    fit_time = lpme_time
  )

  pme_out <- list(
    pme = pme_result_list,
    reconstructions = pme_reconstructions,
    fit_time = pme_time_list
  )

  pc_out <- list(
    pc = pc_result_list,
    reconstructions = pc_reconstructions,
    fit_time = pc_time_list
  )

  model_out <- list(
    lpme = lpme_out,
    pme = pme_out,
    pc = pc_out
  )

  if (partition == TRUE) {
    df_numbered <- data$df_observed |>
      mutate(row_num = row_number())
    mat_part1 <- df_numbered |>
      filter(X1 > 0)
    mat_part2 <- df_numbered |>
      filter(X1 <= 0)
    part1_order <- mat_part1$row_num |>
      unlist()
    part2_order <- mat_part2$row_num |>
      unlist()

    mat_part1 <- mat_part1 |>
      select(time, contains("X")) |>
      as.matrix()
    mat_part2 <- mat_part2 |>
      select(time, contains("X")) |>
      as.matrix()

    lpme_part_results <- replicate(2, NULL, simplify = FALSE)
    pme_part_results <- replicate(
      2,
      replicate(length(time_values), NULL, simplify = FALSE),
      simplify = FALSE
    )
    pc_part_results <- replicate(
      2,
      replicate(length(time_values), NULL, simplify = FALSE),
      simplify = FALSE
    )

    lpme_part_reconstructions <- replicate(2, NULL, simplify = FALSE)
    pme_part_reconstructions <- replicate(
      2,
      replicate(length(time_values), NULL, simplify = FALSE),
      simplify = FALSE
    )
    pc_part_reconstructions <- replicate(
      2,
      replicate(length(time_values), NULL, simplify = FALSE),
      simplify = FALSE
    )

    lpme_part_times <- replicate(2, NULL, simplify = FALSE)
    pme_part_times <- replicate(
      2,
      replicate(length(time_values), NULL, simplify = FALSE),
      simplify = FALSE
    )
    pc_part_times <- replicate(
      2,
      replicate(length(time_values), NULL, simplify = FALSE),
      simplify = FALSE
    )

    tic()
    lpme_part_results[[1]] <- lpme(
      mat_part1,
      d,
      verbose = FALSE,
      print_plots = FALSE
    )
    lpme_part_times[[1]] <- toc()
    lpme_part_reconstructions[[1]] <- calculate_lpme_reconstructions(
      lpme_part_results[[1]],
      mat_part1
    )
    tic()
    lpme_part_results[[2]] <- lpme(
      mat_part2,
      d,
      verbose = FALSE,
      print_plots = FALSE
    )
    lpme_part_times[[2]] <- toc()
    lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
      lpme_part_results[[2]],
      mat_part2
    )

    for (time_idx in seq_along(time_values)) {
      temp_data_part1 <- mat_part1[mat_part1[, 1] == time_values[time_idx], -1]
      temp_data_part2 <- mat_part2[mat_part2[, 1] == time_values[time_idx], -1]

      tic()
      pme_part_results[[1]][[time_idx]] <- pme(temp_data_part1, d = d)
      pme_part_times[[1]][[time_idx]] <- toc()
      pme_part_reconstructions[[1]][[
        time_idx
      ]] <- calculate_pme_reconstructions(
        pme_part_results[[1]][[time_idx]],
        temp_data_part1
      )
      tic()
      pme_part_results[[2]][[time_idx]] <- pme(temp_data_part2, d = d)
      pme_part_times[[2]][[time_idx]] <- toc()
      pme_part_reconstructions[[2]][[
        time_idx
      ]] <- calculate_pme_reconstructions(
        pme_part_results[[2]][[time_idx]],
        temp_data_part2
      )

      pme_part_reconstructions[[1]][[time_idx]] <- cbind(
        time_values[time_idx],
        pme_part_reconstructions[[1]][[time_idx]]
      )
      pme_part_reconstructions[[2]][[time_idx]] <- cbind(
        time_values[time_idx],
        pme_part_reconstructions[[2]][[time_idx]]
      )

      if (d == 1) {
        principal_curves_pt1 <- list()
        pc_times_pt1 <- list()
        pc_error_pt1 <- vector()
        for (smoother_idx in seq_along(smoothing_options)) {
          tic()
          principal_curves_pt1[[smoother_idx]] <- principal_curve(
            temp_data_part1,
            smoother = smoothing_options[[smoother_idx]]
          )
          pc_times_pt1[[smoother_idx]] <- toc()
          pc_error_pt1[smoother_idx] <- principal_curves_pt1[[
            smoother_idx
          ]]$dist
        }
        opt_principal_curve_pt1 <- which.min(pc_error_pt1)
        pc_part_results[[1]][[time_idx]] <- principal_curves_pt1[[
          opt_principal_curve_pt1
        ]]
        pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          pc_part_results[[1]][[time_idx]]$s
        )
        pc_part_times[[1]][[time_idx]] <- pc_times_pt1[[
          opt_principal_curve_pt1
        ]]

        principal_curves_pt2 <- list()
        pc_times_pt2 <- list()
        pc_error_pt2 <- vector()
        for (smoother_idx in seq_along(smoothing_options)) {
          tic()
          principal_curves_pt2[[smoother_idx]] <- principal_curve(
            temp_data_part2,
            smoother = smoothing_options[[smoother_idx]]
          )
          pc_times_pt2[[smoother_idx]] <- toc()
          pc_error_pt2[smoother_idx] <- principal_curves_pt2[[
            smoother_idx
          ]]$dist
        }
        opt_principal_curve_pt2 <- which.min(pc_error_pt2)
        pc_part_results[[2]][[time_idx]] <- principal_curves_pt2[[
          opt_principal_curve_pt2
        ]]
        pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          pc_part_results[[2]][[time_idx]]$s
        )
        pc_part_times[[2]][[time_idx]] <- pc_times_pt2[[
          opt_principal_curve_pt2
        ]]
      } else if (d == 2) {
        tic()
        principal_surface_part1 <- prinSurf(temp_data_part1)
        pc_part_times[[1]][[time_idx]] <- toc()
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        pc_part_results[[1]][[time_idx]] <- principal_surface_part1[[
          opt_surface_part1 + 2
        ]]
        pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )

        tic()
        principal_surface_part2 <- prinSurf(temp_data_part2)
        pc_part_times[[2]][[time_idx]] <- toc()
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        pc_part_results[[2]][[time_idx]] <- principal_surface_part2[[
          opt_surface_part2 + 2
        ]]
        pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
      }
    }

    lpme_part_reconstructions[[1]] <- cbind(
      part1_order,
      lpme_part_reconstructions[[1]]
    )
    lpme_part_reconstructions[[2]] <- cbind(
      part2_order,
      lpme_part_reconstructions[[2]]
    )
    lpme_part_reconstructions <- reduce(lpme_part_reconstructions, rbind)
    lpme_part_reconstructions <- lpme_part_reconstructions[
      order(lpme_part_reconstructions[, 1]),
    ]
    lpme_part_reconstructions <- lpme_part_reconstructions[, -1]

    pme_part_reconstructions <- map(
      seq_along(pme_part_reconstructions),
      ~ reduce(pme_part_reconstructions[[.x]], rbind)
    )
    pme_part_reconstructions[[1]] <- cbind(
      part1_order,
      pme_part_reconstructions[[1]]
    )
    pme_part_reconstructions[[2]] <- cbind(
      part2_order,
      pme_part_reconstructions[[2]]
    )
    pme_part_reconstructions <- reduce(pme_part_reconstructions, rbind)
    pme_part_reconstructions <- pme_part_reconstructions[
      order(pme_part_reconstructions[, 1]),
    ]
    pme_part_reconstructions <- pme_part_reconstructions[, -1]

    pc_part_reconstructions <- map(
      seq_along(pc_part_reconstructions),
      ~ reduce(pc_part_reconstructions[[.x]], rbind)
    )
    pc_part_reconstructions[[1]] <- cbind(
      part1_order,
      pc_part_reconstructions[[1]]
    )
    pc_part_reconstructions[[2]] <- cbind(
      part2_order,
      pc_part_reconstructions[[2]]
    )
    pc_part_reconstructions <- reduce(pc_part_reconstructions, rbind)
    pc_part_reconstructions <- pc_part_reconstructions[
      order(pc_part_reconstructions[, 1]),
    ]
    pc_part_reconstructions <- pc_part_reconstructions[, -1]

    lpme_part_out <- list(
      lpme = lpme_part_results,
      reconstructions = lpme_part_reconstructions,
      fit_time = lpme_part_times
    )
    pme_part_out <- list(
      pme = pme_part_results,
      reconstructions = pme_part_reconstructions,
      fit_time = pme_part_times
    )
    pc_part_out <- list(
      pc = pc_part_results,
      reconstructions = pc_part_reconstructions,
      fit_time = pc_part_times
    )

    part_out <- list(
      lpme_part = lpme_part_out,
      pme_part = pme_part_out,
      pc_part = pc_part_out
    )
  } else {
    part_out <- list(
      lpme_part = NULL,
      pme_part = NULL,
      pc_part = NULL
    )
  }

  model_out <- c(model_out, part_out)

  return(model_out)
}
