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
library(RhpcBLASctl)
library(tidyverse)

handlers(global = TRUE)
handlers("progress")

options(future.globals.maxSize = 32 * 1024^3)
options(renv.config.sandbox.enabled = FALSE)
options(renv.config.auto.snapshot = FALSE)
cores <- availableCores()

blas_set_num_threads(1)
omp_set_num_threads(1)

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
  list.files(here("code/functions/simulations"), full.names = TRUE),
  source
)

case <- 1
if (case == 1) {
  d <- 1
  D <- 2
  template <- "euclidean"
  min_clusters <- 10
  component_type <- "centers"
  N <- 500
} else if (case == 10) {
  d <- 2
  D <- 3
  template <- "sphere"
  min_clusters <- 50
  component_type <- "centers"
  N <- 2000
}


population_time_change <- 0.25
group2_time_change <- c(0, 0.05, 0.1, 0.25, 0.5)
group_time_change_noise <- c(0, 0.025, 0.05, 0.1)
id_time_change_noise <- c(0, 0.025, 0.05, 0.1)
amplitude_noise <- c(0, 0.025, 0.05, 0.1, 0.25)
repetition <- 1:5

sim_param_grid <- expand_grid(
  group2_time_change,
  group_time_change_noise,
  id_time_change_noise,
  amplitude_noise,
  repetition
)


set.seed(100)

plan(multicore, workers = cores)

sim_results <- foreach(sim_idx = seq_len(nrow(sim_param_grid))) %do%
  {
    print(paste0("Starting Simulation ", sim_idx, " of ", nrow(sim_param_grid)))
    rep_num <- sim_param_grid$repetition[sim_idx]

    test_sim <- simulate_hierarchical_data(
      n_groups = 2,
      n_individuals = 100,
      group_probs = c(0.5, 0.5),
      population_time_change = 0.1,
      group_time_change_diff = c(
        -sim_param_grid$group2_time_change[sim_idx],
        sim_param_grid$group2_time_change[sim_idx]
      ),
      group_time_change_noise = sim_param_grid$group_time_change_noise[sim_idx],
      id_time_change_noise = sim_param_grid$id_time_change_noise[sim_idx],
      population_time_trend = "linear",
      group_time_trends = c("linear", "linear"),
      duration = 5,
      interval = 0.5,
      case = case,
      obs_noise = 0.05,
      amplitude_noise = sim_param_grid$amplitude_noise[sim_idx],
      period_noise = 0.1,
      visit_noise = 0.1,
      N = N
    )

    sim_preprocessed <- preprocess_data(test_sim, case = 1, d = d, D = D)

    # LPME Initialization
    init_data <- sim_preprocessed$data |>
      filter(id == "1") |>
      select(time, X1, X2) |>
      as.matrix()

    init_lpme <- lpme(
      init_data,
      d = d,
      template = template,
      init_type = "centers",
      min_clusters = min_clusters,
      print_plots = FALSE,
      verbose = FALSE
    )

    # Data Reduction

    reduced_data <- sim_reduction(
      sim_preprocessed,
      case,
      min_clusters,
      component_type
    )
    sim_reduced <- reduced_data$sim_reduced
    population_reduced <- reduced_data$population_reduced
    groups_reduced <- reduced_data$groups_reduced

    # Parameterization

    sim_x <- sim_preprocessed$data |>
      select(time, contains("X"), -contains("true")) |>
      as.matrix()
    sim_centers <- sim_reduced |>
      select(time_from_bl, contains("X")) |>
      as.matrix()
    sim_weights <- sim_reduced$weight

    init_params <- calculate_lpme_reconstructions(
      init_lpme,
      sim_centers
    )
    init_embeddings <- init_params$reconstructions
    init_params <- init_params$projections

    population_init_params <- calculate_lpme_reconstructions(
      init_lpme,
      as.matrix(select(population_reduced, time_from_bl, contains("X")))
    )$projections

    group_init_params <- foreach(
      group_idx = seq_along(sim_preprocessed$groups)
    ) %do%
      {
        calculate_lpme_reconstructions(
          init_lpme,
          as.matrix(select(
            groups_reduced[[group_idx]],
            time_from_bl,
            contains("X")
          ))
        )$projections
      }

    # Fitting

    if (plot_progress == TRUE) {
      plot_ly(
        x = init_embeddings[, 1],
        y = init_embeddings[, 2],
        z = init_embeddings[, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)
      )
    }

    blas_set_num_threads(cores)

    additive_model_list <- additive_pme(
      reduced_data = sim_reduced,
      centers = list(sim_centers),
      init_params = list(init_params),
      weights = list(sim_weights),
      groups = sim_reduced$group,
      ids = sim_reduced$id,
      scans = sim_reduced$scan,
      times = sim_reduced$time_from_bl,
      partitions = sim_reduced$partition,
      template = template,
      cores = cores,
      lambda = exp(-10:10),
      gamma = exp(-10:10),
      max_iter = 25,
      verbose = FALSE,
      plot_progress = FALSE
    )

    blas_set_num_threads(1)

    additive_model <- additive_model_list$additive_model
    params <- additive_model_list$params
    center_projections <- additive_model_list$center_projections

    # Projections
    print("Computing full projections")

    projection_list <- sim_projections(
      additive_model,
      sim_preprocessed$data,
      sim_reduced,
      params,
      sim_preprocessed$groups,
      sim_preprocessed$ids,
      d = d,
      D = D,
      template = template,
      cores = cores
    )

    source(here("code/functions/simulations/calc_msd.R"))

    additive_full_msd <- calc_msd(
      sim_preprocessed$data,
      projections = projection_list$projections
    )

    population_data <- sim_preprocessed$data_full |>
      filter(id == "Population")

    population_projections <- sim_population_projections(
      additive_model,
      population_data,
      population_reduced,
      list(population_init_params),
      sim_preprocessed$groups,
      sim_preprocessed$ids,
      d = d,
      D = D,
      template = template,
      cores = cores
    )

    population_msd <- calc_msd(
      population_data,
      projections = population_projections$projections
    )

    group_projections <- foreach(
      group_idx = seq_along(sim_preprocessed$groups)
    ) %do%
      {
        group_data <- sim_preprocessed$data_full |>
          filter(id == paste0("Group ", sim_preprocessed$groups[group_idx]))
        group_reduced <- groups_reduced[[group_idx]]

        sim_group_projections(
          additive_model,
          group_data,
          group_reduced,
          list(group_init_params[[group_idx]]),
          group_idx,
          d = d,
          D = D,
          template = template,
          cores = cores
        )
      }

    group_msd <- foreach(group_idx = seq_along(sim_preprocessed$groups)) %do%
      {
        group_data <- sim_preprocessed$data_full |>
          filter(id == paste0("Group ", sim_preprocessed$groups[group_idx]))

        calc_msd(
          group_data,
          projections = group_projections[[group_idx]]$projections
        )
      }

    # as expected, we encounter issues in the convergence of the HDMDE algorithm when
    # applied to the full dataset

    # pop_data <- sim_preprocessed$data |>
    #   select(time, contains("X"), -contains("true")) |>
    #   as.matrix()
    #
    # population_lpme <- lpme(pop_data, d, verbose = FALSE, print_plots = FALSE)
    # population_projections <- calculate_lpme_reconstructions(
    #   population_lpme,
    #   pop_data
    # )
    # population_embeddings <- population_projections$reconstructions
    # population_params <- population_projections$projections
    #
    # pop_lpme_msd <- calc_msd(
    #   sim_preprocessed$data,
    #   population_embeddings
    # )

    # also encountered HDMDE convergence issues when running on group-wide data
    # group_lpme_list <- foreach(
    #   group_idx = seq_along(sim_preprocessed$groups),
    #   .options.future = list(seed = TRUE)
    # ) %dofuture%
    #   {
    #     group_data <- sim_preprocessed$data |>
    #       filter(group == sim_preprocessed$groups[group_idx]) |>
    #       select(time, contains("X"), -contains("true")) |>
    #       as.matrix()
    #
    #     group_lpme <- lpme(group_data, d, verbose = FALSE, print_plots = FALSE)
    #
    #     group_projections <- calculate_lpme_reconstructions(
    #       group_lpme,
    #       group_data
    #     )
    #     group_embeddings <- group_projections$reconstructions
    #     group_params <- group_projections$projections
    #
    #     group_lpme_msd <- calc_msd(
    #       filter(
    #         sim_preprocessed$data,
    #         group == sim_preprocessed$groups[group_idx]
    #       ),
    #       group_embeddings
    #     )
    #
    #     list(
    #       lpme = group_lpme,
    #       msd = group_lpme_msd,
    #       embeddings = group_embeddings,
    #       params = group_params
    #     )
    #   }

    source(here("code/functions/simulations/functional_anova.R"))
    source(here("code/functions/simulations/display_test_results.R"))
    source(here("code/functions/simulations/functional_permutation_test.R"))

    permutation_test_results <- functional_permutation_test(
      additive_model,
      list(sim_centers),
      list(sim_weights),
      params,
      groups = sim_reduced$group,
      ids = sim_reduced$id,
      scans = sim_reduced$scan,
      times = sim_reduced$time_from_bl,
      partitions = sim_reduced$partition,
      id_groups = test_sim$id_groups,
      n_params = 250,
      template = template,
      alpha = 0.05,
      n_permutations = 1000,
      verbose = TRUE
    )

    f_test_any_rejected <- (rowSums(permutation_test_results$f_test$rejected[[
      1
    ]]) >
      0) |>
      mean()
    f_test_all_rejected <- (rowSums(permutation_test_results$f_test$rejected[[
      1
    ]]) ==
      D) |>
      mean()
    f_test_rejected_pct <- mean(permutation_test_results$f_test$rejected[[1]])

    chisq_test_any_rejected <- (rowSums(permutation_test_results$chisq_test$rejected[[
      1
    ]]) >
      0) |>
      mean()
    chisq_test_all_rejected <- (rowSums(permutation_test_results$chisq_test$rejected[[
      1
    ]]) ==
      D) |>
      mean()
    chisq_test_rejected_pct <- mean(permutation_test_results$chisq_test$rejected[[
      1
    ]])

    l2_norm_test_rejected <- permutation_test_results$l2_norm_test$rejected[[1]]

    permutation_test_plots <- display_test_results(
      permutation_test_results,
      sim_preprocessed$groups
    )

    sim_result_out <- list(
      additive_model = additive_model,
      data = test_sim,
      params = projection_list,
      msd = additive_full_msd,
      population_projections = population_projections,
      population_msd = population_msd,
      group_projections = group_projections,
      group_msd = group_msd,
      anova = permutation_test_results,
      anova_plots = permutation_test_plots,

      f_test_results = list(
        any_rejected = f_test_any_rejected,
        all_rejected = f_test_all_rejected,
        rejected_pct = f_test_rejected_pct
      ),
      chisq_test_rejected = list(
        any_rejected = chisq_test_any_rejected,
        all_rejected = chisq_test_all_rejected,
        rejected_pct = chisq_test_rejected_pct
      ),
      l2_norm_test_rejected = l2_norm_test_rejected,
      group_time_change_diff = sim_param_grid$group2_time_change[sim_idx],
      group_time_change_noise = sim_param_grid$group_time_change_noise[sim_idx],
      id_time_change_noise = sim_param_grid$id_time_change_noise[sim_idx],
      amplitude_noise = sim_param_grid$amplitude_noise[sim_idx],
      sim_repetition = rep_num
    )

    filename <- here(
      paste0(
        "output/simulations/case1/sim_case1_group_time_change_diff_",
        str_pad(
          as.character(100 * sim_param_grid$group2_time_change[sim_idx]),
          3,
          pad = "0"
        ),
        "_noise_",
        str_pad(
          as.character(100 * sim_param_grid$group_time_change_noise[sim_idx]),
          3,
          pad = "0"
        ),
        "_id_time_change_noise_",
        str_pad(
          as.character(100 * sim_param_grid$id_time_change_noise[sim_idx]),
          3,
          pad = "0"
        ),
        "_amplitude_noise_",
        str_pad(
          as.character(100 * sim_param_grid$amplitude_noise[sim_idx]),
          3,
          pad = "0"
        ),
        "_rep_",
        as.character(rep_num),
        ".RDS"
      )
    )

    saveRDS(sim_result_out, filename)
    TRUE
  }
