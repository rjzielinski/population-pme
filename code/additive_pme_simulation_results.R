library(dplyr)
library(ggplot2)
library(here)
library(progressr)
library(readr)

handlers(global = TRUE)
handlers("progress")

sim_files <- list.files(
  here("output/simulations/case1"),
  full.names = TRUE
)

n_simulations <- length(sim_files)

sim_index <- vector()

group_time_change_diff <- vector()
group_time_change_noise <- vector()
id_time_change_noise <- vector()

population_msd <- vector()
group_msd1 <- vector()
group_msd2 <- vector()
full_msd <- vector()

f_param_test_any_rejected <- vector()
f_param_test_all_rejected <- vector()
f_param_test_rejected_pct <- vector()

f_permute_test_any_rejected <- vector()
f_permute_test_all_rejected <- vector()
f_permute_test_rejected_pct <- vector()

chisq_param_test_any_rejected <- vector()
chisq_param_test_all_rejected <- vector()
chisq_param_test_rejected_pct <- vector()

chisq_permute_test_any_rejected <- vector()
chisq_permute_test_all_rejected <- vector()
chisq_permute_test_rejected_pct <- vector()

l2_norm_param_test_rejected <- matrix(nrow = n_simulations, ncol = 2)
l2_norm_param_test_p_value <- matrix(nrow = n_simulations, ncol = 2)
l2_norm_permute_test_rejected <- matrix(nrow = n_simulations, ncol = 2)
l2_norm_permute_test_p_value <- matrix(nrow = n_simulations, ncol = 2)

additive_pme_time <- vector()
init_time <- vector()
reduction_time <- vector()
param_time <- vector()
fitting_time <- vector()
projection_time <- vector()
permutation_time <- vector()

with_progress({
  p <- progressor(n_simulations)
  for (sim_idx in seq_along(sim_files)) {
    sim_result <- readRDS(sim_files[sim_idx])

    sim_index[sim_idx] <- sim_idx
    group_time_change_diff[sim_idx] <- sim_result$group_time_change_diff
    group_time_change_noise[sim_idx] <- sim_result$group_time_change_noise
    id_time_change_noise[sim_idx] <- sim_result$id_time_change_noise

    population_msd[sim_idx] <- sim_result$population_msd
    group_msd1[sim_idx] <- sim_result$group_msd[[1]]
    group_msd2[sim_idx] <- sim_result$group_msd[[2]]
    full_msd[sim_idx] <- sim_result$msd

    f_param_test_any_rejected[
      sim_idx
    ] <- sim_result$f_test_rejected$param_any_rejected
    f_param_test_all_rejected[
      sim_idx
    ] <- sim_result$f_test_rejected$param_all_rejected
    f_param_test_rejected_pct[
      sim_idx
    ] <- sim_result$f_test_rejected$param_rejected_pct

    f_permute_test_any_rejected[
      sim_idx
    ] <- sim_result$f_test_rejected$permute_any_rejected
    f_permute_test_all_rejected[
      sim_idx
    ] <- sim_result$f_test_rejected$permute_all_rejected
    f_permute_test_rejected_pct[
      sim_idx
    ] <- sim_result$f_test_rejected$permute_rejected_pct

    chisq_param_test_any_rejected[
      sim_idx
    ] <- sim_result$chisq_test_rejected$param_any_rejected
    chisq_param_test_all_rejected[
      sim_idx
    ] <- sim_result$chisq_test_rejected$param_all_rejected
    chisq_param_test_rejected_pct[
      sim_idx
    ] <- sim_result$chisq_test_rejected$param_rejected_pct

    chisq_permute_test_any_rejected[
      sim_idx
    ] <- sim_result$chisq_test_rejected$permute_any_rejected
    chisq_permute_test_all_rejected[
      sim_idx
    ] <- sim_result$chisq_test_rejected$permute_all_rejected
    chisq_permute_test_rejected_pct[
      sim_idx
    ] <- sim_result$chisq_test_rejected$permute_rejected_pct

    l2_norm_param_test_rejected[
      sim_idx,
    ] <- sim_result$l2_norm_test_rejected$param_test_rejected

    l2_norm_param_test_p_value[
      sim_idx,
    ] <- sim_result$anova$l2_norm_test$param_p_value[[1]]

    l2_norm_permute_test_rejected[
      sim_idx,
    ] <- sim_result$l2_norm_test_rejected$permute_test_rejected
    l2_norm_permute_test_p_value[
      sim_idx,
    ] <- sim_result$anova$l2_norm_test$permute_p_values[[1]][[1]]

    additive_pme_time[sim_idx] <- sim_result$times$init +
      sim_result$times$reduction +
      sim_result$times$parameterization +
      sim_result$times$fitting +
      sim_result$times$projection

    init_time[sim_idx] <- sim_result$times$init
    reduction_time[sim_idx] <- sim_result$times$reduction
    param_time[sim_idx] <- sim_result$times$parameterization
    fitting_time[sim_idx] <- sim_result$times$fitting
    projection_time[sim_idx] <- sim_result$times$projection
    permutation_time[sim_idx] <- sim_result$times$permutation_test

    p()
  }
})

l2_norm_result_df <- cbind(
  l2_norm_param_test_rejected,
  l2_norm_param_test_p_value,
  l2_norm_permute_test_rejected,
  l2_norm_permute_test_p_value
) |>
  as_tibble(.name_repair = "minimal")

names(l2_norm_result_df) <- c(
  "l2_norm_param_test_rejected_x",
  "l2_norm_param_test_rejected_y",
  "l2_norm_param_test_p_value_x",
  "l2_norm_param_test_p_value_y",
  "l2_norm_permute_test_rejected_x",
  "l2_norm_permute_test_rejected_y",
  "l2_norm_permute_test_p_value_x",
  "l2_norm_permute_test_p_value_y"
)

l2_norm_result_df <- l2_norm_result_df |>
  mutate_at(vars(contains("rejected")), as.logical)


sim_result_df <- tibble(
  sim_index,
  group_time_change_diff,
  group_time_change_noise,
  id_time_change_noise,
  population_msd,
  group_msd1,
  group_msd2,
  full_msd,
  f_permute_test_any_rejected,
  f_permute_test_all_rejected,
  f_permute_test_rejected_pct,
  chisq_permute_test_any_rejected,
  chisq_permute_test_all_rejected,
  chisq_permute_test_rejected_pct
) |>
  bind_cols(l2_norm_result_df) |>
  bind_cols(
    tibble(
      additive_pme_time,
      init_time,
      reduction_time,
      param_time,
      fitting_time,
      projection_time,
      permutation_time
    )
  )

write_csv(sim_result_df, here("output/simulation_rejection_results.csv"))
