library(doParallel)
library(doFuture)
library(foreach)
library(here)
library(listenv)
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

# read_data(n_individuals = 100)
read_data()

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
  ssd_prev <- ssd
  additive_model_old <- additive_model
  params_prev <- params

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

projection_list <- final_projections(
  additive_model,
  lhipp_surface,
  lhipp_surface_red,
  params,
  partition_values,
  group_values,
  id_values
)

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
saveRDS(lhipp_test_out, "output/test_lhipp_additive_model.RDS")
