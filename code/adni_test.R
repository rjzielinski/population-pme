library(doParallel)
library(doFuture)
library(foreach)
library(future.mirai)
library(here)
library(listenv)
library(plotly)
library(pme)
library(progressr)
library(Rfast)
library(tidyverse)

handlers(global = TRUE)
handlers("progress")

options(future.globals.maxSize = 32 * 1024^3)
options(renv.config.sandbox.enabled = FALSE)
options(renv.config.auto.snapshot = FALSE)
cores <- detectCores() - 4

plot_progress <- TRUE

source("code/functions/fit_weighted_spline.R")
source("code/functions/print_SSD.R")
source("code/functions/fit_adni_additive_model.R")
source("code/functions/plot_id_pred.R")
source("code/functions/plot_additive_model.R")
source("code/functions/calculate_pme_reconstructions.R")
source("code/functions/calculate_lpme_reconstructions.R")

map(
  list.files(here("code/functions/adni_modeling"), full.names = TRUE),
  source
)

plan(multicore, workers = cores)

ssd_ratio_threshold <- 5
verbose <- TRUE

epsilon <- 0.05
max_iter <- 100

read_data(n_individuals = 100)

# INITIALIZATION

lhipp_init_lpme <- lpme_initialization(
  data = lhipp_surface,
  ids = lhipp_ids,
  d = 2,
  min_clusters = 50,
  print_plots = FALSE,
  verbose = FALSE
)

# DATA REDUCTION

lhipp_reduced <- list()
scan_list <- unique(lhipp_scans)


lhipp_surface_red <- data_reduction(
  lhipp_surface,
  lhipp_groups,
  lhipp_ids,
  lhipp_scans
)

prepare_data(lhipp_surface, lhipp_surface_red)

# PARAMETERIZATION

init_params <- get_init_params(lhipp_init_lpme, lhipp_centers)

partition_values <- unique(lhipp_partition)
partition_values <- partition_values[order(partition_values)]

if (plot_progress == TRUE) {
  init_embeddings <- list()
  for (partition_idx in seq_along(partition_values)) {
    init_embeddings[[partition_idx]] <- future_map(
      seq_len(nrow(init_params[[partition_idx]])),
      ~ lhipp_init_lpme[[partition_idx]]$embedding_map(unlist(init_params[[
        partition_idx
      ]][.x, ])),
      .options = furrr_options(seed = TRUE)
    ) |>
      reduce(rbind)
  }
  init_embeddings_full <- reduce(init_embeddings, rbind)

  plot_ly(
    x = init_embeddings_full[, 2],
    y = init_embeddings_full[, 3],
    z = init_embeddings_full[, 4],
    frame = init_embeddings_full[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  )
}

group_vals <- unique(lhipp_surface_red$group)[order(unique(
  lhipp_surface_red$group
))]

# Fit initial additive model estimate
additive_model <- list()
for (partition_idx in seq_along(partition_values)) {
  additive_model[[partition_idx]] <- fit_adni_additive_model(
    centers = lhipp_centers[[partition_idx]],
    params = init_params[[partition_idx]],
    weights = lhipp_weights[[partition_idx]],
    lambda = exp(-20:5),
    gamma = exp(-20:5),
    groups = filter(
      lhipp_surface_red,
      partition == partition_values[partition_idx]
    )$group,
    ids = filter(
      lhipp_surface_red,
      partition == partition_values[partition_idx]
    )$id,
    scans = filter(
      lhipp_surface_red,
      partition == partition_values[partition_idx]
    )$scan,
    times = filter(
      lhipp_surface_red,
      partition == partition_values[partition_idx]
    )$time_from_bl,
    epsilon = 0.05,
    max_iter = 250,
    cores = cores
  )
}

if (plot_progress == TRUE) {
  plot_additive_model(
    additive_model,
    init_params,
    group_vals,
    margin = 0.02
  )
}

print("Updating parameters")

param_list <- update_params(
  additive_model,
  init_params,
  lhipp_surface_red,
  lhipp_centers,
  partition_values,
  id_values
)

params <- param_list$params
center_projections <- param_list$center_projections

ssd <- calc_ssd(lhipp_centers, center_projections, partition_values)
ssd_ratio <- 10 * epsilon

n <- 1

ssd_vec <- vector()
ssd_vec[n] <- ssd


while (
  (ssd_ratio > epsilon) &
    (ssd_ratio <= ssd_ratio_threshold) &
    (n <= (max_iter - 1))
) {
  ssd_prev <- ssd
  additive_model_old <- additive_model
  params_prev <- params

  print("Fitting additive models")

  for (partition_idx in seq_along(partition_values)) {
    additive_model[[partition_idx]] <- fit_adni_additive_model(
      centers = lhipp_centers[[partition_idx]],
      params = init_params[[partition_idx]],
      weights = lhipp_weights[[partition_idx]],
      lambda = exp(-20:5),
      gamma = exp(-20:5),
      groups = filter(
        lhipp_surface_red,
        partition == partition_values[partition_idx]
      )$group,
      ids = filter(
        lhipp_surface_red,
        partition == partition_values[partition_idx]
      )$id,
      scans = filter(
        lhipp_surface_red,
        partition == partition_values[partition_idx]
      )$scan,
      times = filter(
        lhipp_surface_red,
        partition == partition_values[partition_idx]
      )$time_from_bl,
      epsilon = 0.05,
      max_iter = 250,
      cores = cores
    )
  }

  if (plot_progress == TRUE) {
    fig <- plot_additive_model(
      additive_model,
      params,
      group_vals,
      margin = 0.02
    )
    print(fig)
  }

  print("Updating parameters")

  param_list <- update_params(
    additive_model,
    params,
    lhipp_surface_red,
    lhipp_centers,
    partition_values,
    id_values
  )

  params <- param_list$params
  center_projections <- param_list$center_projections

  ssd <- calc_ssd(lhipp_centers, center_projections, partition_values)
  ssd_ratio <- abs(ssd - ssd_prev) / ssd_prev

  n <- n + 1
  ssd_vec[n] <- ssd

  if (verbose == TRUE) {
    print_ssd(ssd, ssd_ratio, n)
  }
}

if (plot_progress == TRUE) {
  plot_additive_model(
    additive_model,
    params,
    group_vals,
    margin = 0.02
  )
}

# opacity <- 0.75
# line_width <- 2
#
# library(reticulate)
# use_condaenv("lpme")
# pv <- import("pyvista")
#
# time_points <- sort(unique(population_embeddings[, 1]))
#
# palette <- colorRampPalette(hcl.colors(8, palette = "viridis"))
# palette_colors <- palette(151)
# time_colors <- palette_colors[1:101][ceiling((time_points + 1e-10) * 10)]
#
# lhipp_plot <- pv$Plotter(shape = tuple(2L, 2L))
#
# for (time_idx in seq_along(time_points)) {
#   if (time_points[time_idx] %% 1 == 0) {
#     print(time_points[time_idx])
#
#     temp_pop <- population_embeddings[
#       population_embeddings[, 1] == time_points[time_idx],
#       -1
#     ]
#     temp_groups <- map(
#       group_embeddings,
#       ~ .x[.x[, 1] == time_points[time_idx], -1]
#     )
#
#     pop_projection3 <- mesh_projection(temp_pop, axis = 3)
#     group_projections <- map(temp_groups, ~ mesh_projection(.x, axis = 3))
#
#     lhipp_plot$subplot(0L, 0L)
#     lhipp_plot$add_mesh(
#       pop_projection3$boundary,
#       show_edges = TRUE,
#       color = time_colors[time_idx],
#       line_width = line_width,
#       opacity = opacity,
#       label = paste0("Time = ", round(time_points[time_idx], 2))
#     )
#
#     if (time_idx == 1) {
#       lhipp_plot$add_text(
#         "Population",
#         position = "upper_left",
#         font_size = 12
#       )
#     }
#
#     lhipp_plot$subplot(0L, 1L)
#     lhipp_plot$add_mesh(
#       group_projections[[1]]$boundary,
#       show_edges = TRUE,
#       color = time_colors[time_idx],
#       line_width = line_width,
#       opacity = opacity,
#       label = paste0("Time = ", round(time_points[time_idx], 2))
#     )
#
#     if (time_idx == 1) {
#       lhipp_plot$add_text(
#         group_vals[1],
#         position = "upper_left",
#         font_size = 12
#       )
#     }
#
#     lhipp_plot$subplot(1L, 0L)
#     lhipp_plot$add_mesh(
#       group_projections[[2]]$boundary,
#       show_edges = TRUE,
#       color = time_colors[time_idx],
#       line_width = line_width,
#       opacity = opacity,
#       label = paste0("Time = ", round(time_points[time_idx], 2))
#     )
#
#     if (time_idx == 1) {
#       lhipp_plot$add_text(
#         group_vals[2],
#         position = "upper_left",
#         font_size = 12
#       )
#     }
#
#     lhipp_plot$subplot(1L, 1L)
#     lhipp_plot$add_mesh(
#       group_projections[[3]]$boundary,
#       show_edges = TRUE,
#       color = time_colors[time_idx],
#       line_width = line_width,
#       opacity = opacity,
#       label = paste0("Time = ", round(time_points[time_idx], 2))
#     )
#
#     if (time_idx == 1) {
#       lhipp_plot$add_text(
#         group_vals[3],
#         position = "upper_left",
#         font_size = 12
#       )
#     }
#
#     lhipp_plot$link_views()
#     lhipp_plot$camera_position <- "xy"
#   } else {
#     next
#   }
# }
#
# lhipp_plot$show()
#
# lhipp_plot$screenshot(
#   "output/additive_mod_silhouette_plot.png"
# )
#

# calculate final MSD and projections?

print("Computing full projections")
projection_list <- final_projections(
  additive_model,
  lhipp_surface,
  lhipp_surface_red,
  params,
  partition_values,
  group_values,
  id_values,
  d = ncol(params[[1]]) - 1,
  D = ncol(lhipp_centers[[1]]) - 1,
  cores = cores
)

print("Calculating MSD")
msd <- calc_msd(lhipp_surface, projections = projection_list$projections)

lhipp_test_out <- list(
  model = additive_model,
  data = lhipp_surface,
  reduced_data = lhipp_surface_red,
  params = params,
  projections = projection_list,
  group_values = group_values,
  id_values = id_values,
  partition_values = partition_values,
  msd = msd
)
saveRDS(lhipp_test_out, "output/lhipp_additive_model_100.RDS")

# lhipp_test_out <- readRDS("output/test_lhipp_additive_model.RDS")

# lhipp_test_projections <- final_projections(
#   lhipp_test_out$model,
#   surface_data = lhipp_test_out$data,
#   reduced_data = lhipp_test_out$reduced_data,
#   params = lhipp_test_out$params,
#   partition_values = lhipp_test_out$partition_values,
#   group_values = lhipp_test_out$group_values,
#   id_values = lhipp_test_out$id_values
# )

#
# data <- lhipp_test_out$data
#
# groups <- lhipp_test_out$group_values
# ids <- lhipp_test_out$id_values
#
# id_groups <- map(
#   ids,
#   ~ filter(data, subid == .x) |>
#     pull(Group) |>
#     unique()
# ) |>
#   reduce(c)
#
# plan(multicore, workers = cores)
#
# f_test_results <- functional_anova(
#   lhipp_test_out$model,
#   lhipp_test_out$params,
#   groups,
#   ids,
#   id_groups,
#   n_params = 1000,
#   test_type = "f_type",
#   alpha = 0.05
# )
#
# f_test_bootstrap_results <- functional_anova(
#   lhipp_test_out$model,
#   lhipp_test_out$params,
#   groups,
#   ids,
#   id_groups,
#   n_params = 1000,
#   test_type = "f_type",
#   alpha = 0.05,
#   bootstrap = TRUE,
#   n_bootstrap = 1000
# )
#
# chisq_test_results <- functional_anova(
#   lhipp_test_out$model,
#   lhipp_test_out$params,
#   groups,
#   ids,
#   id_groups,
#   n_params = 1000,
#   test_type = "chisq_type",
#   alpha = 0.05
# )
#
# chisq_test_bootstrap_results <- functional_anova(
#   lhipp_test_out$model,
#   lhipp_test_out$params,
#   groups,
#   ids,
#   id_groups,
#   n_params = 1000,
#   test_type = "chisq_type",
#   alpha = 0.05,
#   bootstrap = TRUE,
#   n_bootstrap = 1000
# )
#
# l2_norm_results <- functional_anova(
#   lhipp_test_out$model,
#   lhipp_test_out$params,
#   groups,
#   ids,
#   id_groups,
#   n_params = 1000,
#   test_type = "l2_norm",
#   alpha = 0.05
# )
#
# l2_norm_bootstrap_results <- functional_anova(
#   lhipp_test_out$model,
#   lhipp_test_out$params,
#   groups,
#   ids,
#   id_groups,
#   n_params = 1000,
#   test_type = "l2_norm",
#   alpha = 0.05,
#   bootstrap = TRUE,
#   n_bootstrap = 1000
# )
#
# display_test_results(pointwise_bootstrap_results, groups)
#
# display_test_results(f_test_results, groups)
