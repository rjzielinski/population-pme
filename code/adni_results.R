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

lhipp_out200 <- readRDS("output/lhipp_additive_model_200.RDS")

silhouette_plot(
  additive_model = lhipp_out200$model,
  surface_data = lhipp_out200$data,
  projections = lhipp_out200$projections,
  axis = 1,
  filename = here("output/lhipp200_silhouette_plot_ax1.png")
)

silhouette_plot(
  additive_model = lhipp_out200$model,
  surface_data = lhipp_out200$data,
  projections = lhipp_out200$projections,
  axis = 2,
  filename = here("output/lhipp200_silhouette_plot_ax2.png")
)

silhouette_plot(
  additive_model = lhipp_out200$model,
  surface_data = lhipp_out200$data,
  projections = lhipp_out200$projections,
  axis = 3,
  filename = here("output/lhipp200_silhouette_plot_ax3.png")
)

data <- lhipp_out200$data

groups <- lhipp_out200$group_values
ids <- lhipp_out200$id_values

id_groups <- map(
  ids,
  ~ filter(data, subid == .x) |>
    pull(Group) |>
    unique()
) |>
  reduce(c)

# plan(multicore, workers = cores / 2)
plan(sequential)

f_test_results <- functional_anova(
  lhipp_out200$model,
  lhipp_out200$params,
  groups,
  ids,
  id_groups,
  n_params = 500,
  test_type = "f_type",
  alpha = 0.05
)
display_test_results(f_test_results, groups)

f_test_bootstrap_results <- functional_anova(
  lhipp_out200$model,
  lhipp_out200$params,
  groups,
  ids,
  id_groups,
  n_params = 500,
  test_type = "f_type",
  alpha = 0.05,
  bootstrap = TRUE,
  n_bootstrap = 1000
)

chisq_test_results <- functional_anova(
  lhipp_out200$model,
  lhipp_out200$params,
  groups,
  ids,
  id_groups,
  n_params = 1000,
  test_type = "chisq_type",
  alpha = 0.05
)

chisq_test_bootstrap_results <- functional_anova(
  lhipp_out200$model,
  lhipp_out200$params,
  groups,
  ids,
  id_groups,
  n_params = 1000,
  test_type = "chisq_type",
  alpha = 0.05,
  bootstrap = TRUE,
  n_bootstrap = 1000
)

l2_norm_results <- functional_anova(
  lhipp_out200$model,
  lhipp_out200$params,
  groups,
  ids,
  id_groups,
  n_params = 1000,
  test_type = "l2_norm",
  alpha = 0.05
)

l2_norm_bootstrap_results <- functional_anova(
  lhipp_out200$model,
  lhipp_out200$params,
  groups,
  ids,
  id_groups,
  n_params = 1000,
  test_type = "l2_norm",
  alpha = 0.05,
  bootstrap = TRUE,
  n_bootstrap = 1000
)

display_test_results(pointwise_bootstrap_results, groups)
