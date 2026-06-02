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

source("code/functions/additive_pme.R")
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


ssd_ratio_threshold <- 5
verbose <- TRUE

epsilon <- 0.05
max_iter <- 100

read_data(
  n_partitions = 1,
  ad_cn_ratio = 1,
  ad_mci_ratio = 1
)

# INITIALIZATION

lhipp_init_lpme <- lpme_initialization(
  data = lhipp_surface,
  ids = lhipp_ids,
  n_init = 1,
  d = 2,
  template = "sphere",
  init_type = "centers",
  min_clusters = 75,
  print_plots = FALSE,
  verbose = FALSE
)

# DATA REDUCTION

plan(multicore, workers = cores)

lhipp_reduced <- list()
scan_list <- unique(lhipp_scans)


lhipp_surface_red <- data_reduction(
  lhipp_surface,
  lhipp_groups,
  lhipp_ids,
  lhipp_scans,
  min_clusters = 75,
  component_type = "centers"
)

prepare_data(lhipp_surface, lhipp_surface_red)

# PARAMETERIZATION

init_params <- get_init_params(lhipp_init_lpme, lhipp_centers)
init_embeddings <- init_params$reconstructions
init_params <- init_params$params

partition_values <- unique(lhipp_partition)
partition_values <- partition_values[order(partition_values)]

group_vals <- sort(unique(lhipp_surface_red$group))

if (plot_progress == TRUE) {
  plot_ly(
    x = init_embeddings[, 2],
    y = init_embeddings[, 3],
    z = init_embeddings[, 4],
    frame = init_embeddings[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  )
}

additive_model_list <- additive_pme(
  reduced_data = lhipp_surface_red,
  centers = lhipp_centers,
  init_params = init_params,
  weights = lhipp_weights,
  groups = lhipp_surface_red$group,
  ids = lhipp_surface_red$id,
  scans = lhipp_surface_red$scan,
  times = lhipp_surface_red$time_from_bl,
  partitions = lhipp_surface_red$partition,
  template = "sphere",
  cores = cores,
  plot_progress = TRUE,
)

additive_model <- additive_model_list$additive_model
params <- additive_model_list$params
center_projections <- additive_model_list$center_projections

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
  template = "sphere",
  cores = cores
)

print("Calculating MSD")
msd <- calc_msd(lhipp_surface, projections = projection_list$projections)

lhipp_test_out <- list(
  model = additive_model,
  data = lhipp_surface,
  reduced_data = lhipp_surface_red,
  projections = projection_list,
  params = params,
  center_projections = center_projections,
  group_values = group_values,
  id_values = id_values,
  partition_values = partition_values,
  msd = msd,
  template = "sphere"
)
saveRDS(lhipp_test_out, "output/lhipp_additive_matched_111.RDS")
