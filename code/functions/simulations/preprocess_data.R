preprocess_data <- function(sim_data, case, d, D) {
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(pracma, quietly = TRUE, warn.conflicts = FALSE)

  df <- sim_data$data

  ##### PREPROCESSING *****

  if (case == 1) {
    df_centers <- df |>
      group_by(id) |>
      summarize(
        mean_x1 = mean(X1),
        mean_x2 = mean(X2),
        mean_x_true1 = mean(X_true1),
        mean_x_true2 = mean(X_true2),
        max_x1 = max(abs(X1 - mean_x1)),
        max_x2 = max(abs(X2 - mean_x2))
      )

    df_bl <- df |>
      group_by(id) |>
      arrange(time) |>
      summarize(duration = max(time))

    df <- df |>
      full_join(df_centers, by = c("id")) |>
      full_join(df_bl, by = "id") |>
      mutate(
        X1 = (X1 - mean_x1) / max_x1,
        X2 = (X2 - mean_x2) / max_x2,
        X_true1 = (X_true1 - mean_x_true1) / max_x1,
        X_true2 = (X_true2 - mean_x_true2) / max_x2
      ) |>
      select(-contains("mean_"), -contains("max_"), -duration)
  } else if (case == 10) {
    df_centers <- df |>
      group_by(id, scan)
  }

  df_full <- df
  df <- df |>
    filter(id != "Population") |>
    filter(!(grepl("Group", id)))

  out_list <- list(
    data = df,
    data_full = df_full,
    groups = sim_data$groups,
    ids = sim_data$ids,
    scans = sim_data$scans
  )
}
