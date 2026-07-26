library(doFuture, quietly = TRUE)
library(foreach, quietly = TRUE)
library(furrr, quietly = TRUE)
library(future, quietly = TRUE)
library(here, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(pme, quietly = TRUE)
library(progressr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(RhpcBLASctl, quietly = TRUE)

source(here("code/functions/adni_modeling/functional_permutation_test.R"))
args <- commandArgs(trailingOnly = TRUE)
permute_idx <- as.integer(args[1])

threads <- 4
RhpcBLASctl::blas_set_num_threads(threads)
RhpcBLASctl::omp_set_num_threads(threads)

n_params <- 1000
template <- "sphere"
mode <- "additive_embeddings"

lambda <- exp(-15:5)
gamma <- exp(-15:5)

model_out <- readRDS("output/lhipp_additive_model_matched_111_sub90.RDS")

additive_model <- model_out$additive_model
centers <- model_out$centers
weights <- model_out$weights
params <- model_out$params
groups <- model_out$groups
ids <- model_out$ids
scans <- model_out$scans
times <- model_out$times
partitions <- model_out$partitions
id_groups <- model_out$id_groups

n_partitions <- length(additive_model)
n_groups <- length(additive_model[[1]]$group_embeddings)
n_individuals <- length(additive_model[[1]]$id_embeddings)

group_vals <- sort(unique(groups))
id_vals <- sort(unique(ids))
scan_vals <- sort(unique(scans))
partition_vals <- sort(unique(partitions))

param_grids <- calc_param_grids(
  params,
  n_params,
  interval = 0.25,
  template = template
)

param_interval_vols <- vector()

for (partition_idx in seq_len(n_partitions)) {
  param_interval_length <- vector()

  if (template == "euclidean") {
    for (dim_idx in seq_len(ncol(param_grids[[partition_idx]]))) {
      values <- param_grids[[partition_idx]][, dim_idx] |>
        unique() |>
        sort()

      param_interval_length[dim_idx] <- values[2] - values[1]
    }

    param_interval_vols[partition_idx] <- prod(param_interval_length)
  } else if (template == "sphere") {
    param_interval_vols[partition_idx] <- (2 * pi) * pi / n_params
  }
}

set.seed(100 + permute_idx)

RhpcBLASctl::blas_set_num_threads(threads)
RhpcBLASctl::omp_set_num_threads(threads)

permute_id_groups <- id_groups |>
  sample(size = n_individuals)

permute_groups <- vector(length = length(groups))
for (id_idx in seq_along(id_vals)) {
  id_set <- ids == id_vals[id_idx]
  id_group_val <- permute_id_groups[id_idx]
  permute_groups[id_set] <- id_group_val
}

permute_model <- list()
for (partition_idx in seq_len(n_partitions)) {
  partition_groups <- permute_groups[
    partitions == partition_vals[partition_idx]
  ]
  partition_ids <- ids[partitions == partition_vals[partition_idx]]
  partition_scans <- scans[partitions == partition_vals[partition_idx]]
  partition_times <- times[partitions == partition_vals[partition_idx]]

  permute_model[[partition_idx]] <- fit_permutation_models(
    additive_model[[partition_idx]]$population_embedding,
    centers[[partition_idx]],
    params[[partition_idx]],
    weights[[partition_idx]],
    additive_model[[partition_idx]]$param_grid,
    lambda = lambda,
    gamma = gamma,
    groups = partition_groups,
    ids = partition_ids,
    scans = partition_scans,
    times = partition_times,
    template = template,
    epsilon = 0.025,
    max_iter = 25,
    cores = 1,
    plot_progress = FALSE
  )
}

if (mode == "additive_embeddings") {
  permute_embeddings <- calc_embeddings(
    permute_model,
    param_grids,
    n_groups,
    n_individuals,
    permute_id_groups
  )

  gc()

  permute_sum_squares <- calc_sum_squares(
    permute_embeddings,
    param_grids,
    permute_groups,
    ids
  )
} else if (mode == "sample_mean") {
  permute_embeddings <- calc_embeddings_sample_means(
    permute_model,
    param_grids,
    n_groups,
    n_individuals,
    permute_id_groups
  )

  gc()

  permute_sum_squares <- calc_sum_squares_sample_means(
    permute_embeddings,
    param_grids,
    permute_groups,
    ids
  )
}

permute_f_stat <- list()
permute_l2_norm_stat <- list()

for (partition_idx in seq_len(n_partitions)) {
  permute_f_stat[[partition_idx]] <- (permute_sum_squares$ssh[[
    partition_idx
  ]] /
    (n_groups - 1)) /
    (permute_sum_squares$sse[[partition_idx]] /
      (n_individuals - n_groups))

  permute_l2_norm_stat[[partition_idx]] <- colSums(
    permute_sum_squares$ssh[[partition_idx]] *
      param_interval_vols[partition_idx]
  )
}

permutation_stats <- list(
  f_test_stat = permute_f_stat,
  l2_norm_test_stat = permute_l2_norm_stat
)

out_dir <- here("output/adni_permutations/")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

saveRDS(
  permutation_stats,
  file = here(paste0(
    "output/adni_permutations/lhipp_permutation_stats_",
    permute_idx,
    ".RDS"
  ))
)
