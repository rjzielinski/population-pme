library(doParallel)
library(doFuture)
library(foreach)
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
cores <- detectCores() - 2

source("code/functions/fit_weighted_spline.R")
source("code/functions/print_SSD.R")
source("code/functions/fit_adni_additive_model.R")
source("code/functions/plot_id_pred.R")
source("code/functions/calculate_pme_reconstructions.R")
source("code/functions/calculate_lpme_reconstructions.R")

map(
  list.files(here("code/functions/adni_modeling"), full.names = TRUE),
  source
)

# plan(multicore, workers = cores)
plan(sequential)

ssd_ratio_threshold <- 5
verbose <- TRUE

epsilon <- 0.05
max_iter <- 100

read_data(n_individuals = 50)

gc()

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

# Fit initial additive model estimate
init_additive_mod <- list()
for (partition_idx in seq_along(partition_values)) {
  init_additive_mod[[partition_idx]] <- fit_adni_additive_model(
    centers = lhipp_centers[[partition_idx]],
    params = init_params[[partition_idx]],
    weights = lhipp_weights[[partition_idx]],
    lambda = exp(-20:5),
    gamma = exp(-20:5),
    k = 5,
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

plan(sequential)

print("Updating parameters")

param_list <- update_params(
  init_additive_mod,
  init_params,
  lhipp_surface_red,
  lhipp_centers,
  partition_values,
  id_values
)

params <- param_list$params
center_projections <- param_list$center_projections

ssd <- calc_ssd(lhipp_centers, center_projections, partition_values)

# embeddings <- init_additive_mod$embeddings

n <- 1
ssd_ratio <- 10 * epsilon

additive_model <- init_additive_mod

ssd_vec <- vector()
ssd_vec[0] <- ssd

while (
  (ssd_ratio > epsilon) &
    (ssd_ratio <= ssd_ratio_threshold) &
    (n <= (max_iter - 1))
) {
  plan(multicore, workers = cores)

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
      k = 5,
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

  plan(sequential)

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

  ssd_vec[n] <- ssd

  n <- n + 1
  if (verbose == TRUE) {
    print_ssd(ssd, ssd_ratio, n)
  }
}

# calculate final MSD and projections?

plan(sequential)

print("Computing full projections")
projection_list <- final_projections(
  additive_model,
  lhipp_surface,
  lhipp_surface_red,
  params,
  partition_values,
  group_values,
  id_values
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
  msd = msd
)
saveRDS(lhipp_test_out, "output/lhipp_additive_model_250.RDS")

lhipp_test_out <- readRDS("output/test_lhipp_additive_model.RDS")

data <- lhipp_test_out$data

groups <- lhipp_test_out$group_values
ids <- lhipp_test_out$id_values
id_groups <- map(
  ids,
  ~ filter(data, subid == .x) |>
    pull(Group) |>
    unique()
) |>
  reduce(c)

f_test_results <- functional_anova_pointwise(
  lhipp_test_out$model,
  lhipp_test_out$params,
  groups,
  ids,
  id_groups,
  n_params = 1000,
  test_type = "f_type",
  alpha = 0.05
)

f_test_group_embeddings <- list()

for (partition_idx in seq_along(lhipp_test_out$model)) {
  partition_group_embeddings <- list()
  partition_params <- lhipp_test_out$params[[partition_idx]]

  for (group_idx in seq_along(groups)) {
    # partition_group_embeddings[[group_idx]] <- map(

    test_embeddings <- map(
      seq_len(nrow(partition_params)),
      ~ {
        cbind(
          partition_params[, 1],
          lhipp_test_out$model[[
            partition_idx
          ]]$population_embedding$embedding_map(unlist(partition_params[
            .x,
          ]))[-1] +
            lhipp_test_out$model[[partition_idx]]$group_embeddings[[
              group_idx
            ]]$embedding_map(unlist(partition_params[.x, ]))[-1]
        )
      }
    ) |>
      reduce(rbind)
  }
}

chisq_test_results <- functional_anova_pointwise(
  lhipp_test_out$model,
  lhipp_test_out$params,
  groups,
  ids,
  id_groups,
  n_params = 1000,
  test_type = "chisq_type",
  alpha = 0.05
)
