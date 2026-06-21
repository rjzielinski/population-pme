simulate_hierarchical_data <- function(
  n_groups,
  n_individuals,
  group_probs,
  population_time_change,
  group_time_change_diff,
  group_time_change_noise,
  id_time_change_noise,
  population_time_trend,
  group_time_trends,
  duration,
  interval,
  case,
  obs_noise,
  amplitude_noise,
  period_noise,
  visit_noise = 0.1,
  N = 1000
) {
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(progressr, quietly = TRUE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  require(tibble, quietly = TRUE, warn.conflicts = FALSE)

  source("code/functions/simulations/simulate_data.R")

  if (is.null(group_probs)) {
    group_probs <- rep(1 / n_groups, n_groups)
  }

  id_groups <- rmultinom(n_individuals, size = 1, prob = group_probs) |>
    apply(MARGIN = 2, FUN = which.max)

  if (is.null(group_time_change_diff)) {
    group_time_change_diff <- rep(0, n_groups)
  }

  group_time_changes <- rnorm(
    n_groups,
    population_time_change + group_time_change_diff,
    group_time_change_noise
  )

  id_time_changes <- rnorm(
    n_individuals,
    group_time_changes[id_groups],
    id_time_change_noise
  )

  if (is.null(group_time_trends)) {
    group_time_trends <- rep("linear", n_groups)
  }

  n_scans <- seq(0, duration, interval) |>
    length()

  scan_counters <- n_scans * (1:n_individuals - 1) + 1

  with_progress({
    p <- progressor(n_individuals)
    simulated_datasets <- foreach(
      id_idx = seq_len(n_individuals),
      .options.future = list(seed = TRUE)
    ) %dofuture%
      {
        id_data <- simulate_data(
          duration,
          interval,
          case,
          obs_noise,
          amplitude_noise,
          period_noise,
          group_time_trends[id_groups[id_idx]],
          id_time_changes[id_idx],
          visit_noise,
          N
        )

        id_df <- id_data$df
        id_times <- unique(id_df$time)

        scan_ids <- scan_counters[id_idx]:(scan_counters[id_idx] +
          length(id_times) -
          1)
        scan_vals <- map(
          id_df$time,
          ~ which(.x == id_times)
        ) |>
          reduce(c)

        id_df$id <- id_idx
        id_df$group <- id_groups[id_idx]
        id_df$scan <- scan_ids[scan_vals]

        p()
        id_df
      }
  })

  sim_df <- do.call(bind_rows, simulated_datasets)

  scan_counter <- max(scan_counters) + n_scans

  population_manifold <- simulate_data(
    duration,
    interval,
    case,
    obs_noise = 0,
    amplitude_noise = 0,
    period_noise = 0,
    population_time_trend,
    population_time_change,
    visit_noise = 0,
    N
  )

  pop_times <- unique(population_manifold$df$time)

  scan_ids <- scan_counter:(scan_counter + length(pop_times) - 1)
  scan_vals <- map(
    population_manifold$df$time,
    ~ which(.x == pop_times)
  ) |>
    reduce(c)

  population_manifold$df <- population_manifold$df |>
    mutate(
      id = "Population",
      group = NA
    )
  population_manifold$df$scan <- scan_ids[scan_vals]

  scan_counter <- scan_counter + length(pop_times)

  true_df <- population_manifold$df

  group_manifolds <- list()
  for (group_idx in seq_len(n_groups)) {
    group_manifolds[[group_idx]] <- simulate_data(
      duration,
      interval,
      case,
      obs_noise = 0,
      amplitude_noise = 0,
      period_noise = 0,
      group_time_trends[group_idx],
      group_time_changes[group_idx],
      visit_noise = 0,
      N
    )

    group_times <- unique(group_manifolds[[group_idx]]$df$time)

    scan_ids <- scan_counter:(scan_counter + length(group_times) - 1)
    scan_vals <- map(
      group_manifolds[[group_idx]]$df$time,
      ~ which(.x == group_times)
    ) |>
      reduce(c)

    group_manifolds[[group_idx]]$df <- group_manifolds[[group_idx]]$df |>
      mutate(
        id = paste0("Group ", group_idx),
        group = group_idx
      )
    group_manifolds[[group_idx]]$df$scan <- scan_ids[scan_vals]

    scan_counter <- scan_counter + length(group_times)

    true_df <- bind_rows(true_df, group_manifolds[[group_idx]]$df)
  }

  groups <- sort(unique(sim_df$group))
  ids <- sort(unique(sim_df$id))
  scans <- sort(unique(sim_df$scan))

  sim_df <- sim_df |>
    mutate(id = as.character(id)) |>
    bind_rows(true_df)

  sim_out <- list(
    data = sim_df,
    population_manifold = population_manifold,
    group_manifolds = group_manifolds,
    groups = groups,
    ids = ids,
    scans = scans,
    id_groups = id_groups
  )

  sim_out
}
