read_data <- function(
  n_individuals = NULL,
  min_visits = 3,
  min_duration = 2,
  n_partitions = 1,
  x_bound = 0,
  y_bound = 0,
  z_bound = 0,
  x_slope = 1,
  scan_exclude_threshold = 0.0025,
  ad_cn_ratio = 1,
  ad_mci_ratio = 1
) {
  # Read and prepare data for anlaysis

  require(cobalt, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(lubridate, quietly = TRUE)
  require(MatchIt, quietly = TRUE)
  require(readr, quietly = TRUE)

  lhipp_surface <- read_csv(here("data/lhipp_surface_fsl.csv"))
  adni_info <- read_csv(here("data/adni_info_full.csv"))

  adni_info <- adni_info |>
    distinct(.keep_all = TRUE)

  lhipp_surface <- lhipp_surface |>
    left_join(adni_info, by = join_by(scan_id == image_id)) |>
    mutate(
      date = decimal_date(date),
      Group = gsub("EMCI", "MCI", Group),
      Group = gsub("LMCI", "MCI", Group),
    )

  adni_n <- lhipp_surface |>
    group_by(subid, scan_id) |>
    tally() |>
    ungroup()
  adni_threshold <- quantile(adni_n$n, scan_exclude_threshold)

  exclude_scans <- adni_n |>
    filter(n < adni_threshold) |>
    select(scan_id) |>
    unlist()

  lhipp_surface <- lhipp_surface |>
    filter(!(scan_id %in% exclude_scans))

  lhipp_surface_centers <- lhipp_surface |>
    group_by(subid, date, scan_id) |>
    summarize(
      mean_x = mean(x),
      mean_y = mean(y),
      mean_z = mean(z),
      max_x = max(abs(x - mean_x)),
      max_y = max(abs(y - mean_y)),
      max_z = max(abs(z - mean_z))
    )

  lhipp_bl <- lhipp_surface |>
    group_by(subid) |>
    arrange(date) |>
    summarize(
      date_bl = first(date),
      max_date = max(date)
    ) |>
    mutate(duration = max_date - date_bl)

  lhipp_surface <- lhipp_surface |>
    full_join(lhipp_surface_centers, by = c("subid", "date", "scan_id")) |>
    full_join(lhipp_bl, by = "subid") |>
    mutate(
      x = (x - mean_x) / max_x,
      y = (y - mean_y) / max_y,
      z = (z - mean_z) / max_z,
      time_from_bl = date - date_bl
    ) |>
    rename(image_id = scan_id)

  if (n_partitions == 2) {
    lhipp_surface <- lhipp_surface |>
      mutate(
        partition = ifelse(z > z_bound, 1, 2)
      )
  } else if (n_partitions == 4) {
    lhipp_surface <- lhipp_surface |>
      mutate(
        partition = case_when(
          y > (-(x_slope * x) + y_bound) & z > z_bound ~ 1,
          y > (-(x_slope * x) + y_bound) & z <= z_bound ~ 2,
          y <= (-(x_slope * x) + y_bound) & z > z_bound ~ 3,
          y <= (-(x_slope * x) + y_bound) & z <= z_bound ~ 4
        )
      )
  } else if (n_partitions == 8) {
    lhipp_surface <- lhipp_surface |>
      mutate(
        partition = case_when(
          x > x_bound & y > y_bound & z > z_bound ~ 1,
          x > x_bound & y > y_bound & z <= z_bound ~ 2,
          x > x_bound & y <= y_bound & z > z_bound ~ 3,
          x > x_bound & y <= y_bound & z <= z_bound ~ 4,
          x <= x_bound & y > y_bound & z > z_bound ~ 5,
          x <= x_bound & y > y_bound & z <= z_bound ~ 6,
          x <= x_bound & y <= y_bound & z > z_bound ~ 7,
          x <= x_bound & y <= y_bound & z <= z_bound ~ 8
        )
      )
  } else {
    lhipp_surface <- lhipp_surface |>
      mutate(partition = 1)
  }

  n_visits <- lhipp_surface |>
    group_by(subid, time_from_bl) |>
    tally() |>
    ungroup() |>
    group_by(subid) |>
    tally() |>
    arrange(n)

  lhipp_surface <- lhipp_surface |>
    filter(subid %in% filter(n_visits, n >= min_visits)$subid) |>
    filter(duration > min_duration)

  id_vec <- unique(lhipp_surface$subid)
  if (!is.null(n_individuals)) {
    n_individuals <- min(n_individuals, length(id_vec))
  } else {
    n_individuals <- length(id_vec)
  }

  include_ids <- id_vec[1:n_individuals]

  lhipp_surface_info <- lhipp_surface |>
    select(subid, Group, Sex, Age, duration) |>
    group_by(subid) |>
    summarize(
      Group = first(Group),
      Sex = first(Sex),
      Age = first(Age),
      Duration = first(duration)
    ) |>
    ungroup() |>
    mutate(Group = as.factor(Group))

  lhipp_surface_info_ad_cn <- lhipp_surface_info |>
    filter(Group %in% c("AD", "CN")) |>
    mutate(
      Group = relevel(Group, ref = "CN"),
      Group = as.numeric(Group) - 1,
      Group = as.logical(Group),
      Sex = as.factor(Sex),
      Sex = as.numeric(Sex)
    )

  lhipp_surface_info_ad_mci <- lhipp_surface_info |>
    filter(Group %in% c("AD", "MCI")) |>
    mutate(
      Group = relevel(Group, ref = "MCI"),
      Group = as.numeric(Group) - 1,
      Group = as.logical(Group),
      Sex = as.factor(Sex),
      Sex = as.numeric(Sex)
    )

  match_ad_cn <- matchit(
    Group ~ Sex + Age + Duration,
    data = lhipp_surface_info_ad_cn,
    method = "nearest",
    distance = "glm",
    ratio = ad_cn_ratio
  )

  match_ad_mci <- matchit(
    Group ~ Sex + Age + Duration,
    data = lhipp_surface_info_ad_mci,
    method = "nearest",
    distance = "glm",
    ratio = ad_mci_ratio
  )

  lhipp_surface_ad_cn_include <- lhipp_surface_info_ad_cn[
    as.logical(match_ad_cn$weights),
  ]

  lhipp_surface_ad_mci_include <- lhipp_surface_info_ad_mci[
    as.logical(match_ad_mci$weights),
  ]

  include_ids <- c(
    lhipp_surface_ad_cn_include$subid,
    lhipp_surface_ad_mci_include$subid
  ) |>
    unique()

  lhipp_surface <- lhipp_surface |>
    filter(subid %in% include_ids)

  group_ids <- lhipp_surface |>
    group_by(Group, subid) |>
    tally() |>
    ungroup() |>
    group_by(Group) |>
    tally()
  include_groups <- group_ids |>
    filter(n >= 5) |>
    select(Group) |>
    reduce(c)

  lhipp_surface <- lhipp_surface |>
    filter(Group %in% include_groups)

  lhipp_groups <- unique(lhipp_surface$Group)
  lhipp_ids <- unique(lhipp_surface$subid)
  lhipp_scans <- unique(lhipp_surface$image_id)
  lhipp_partition <- lhipp_surface$partition

  out_list <- list(
    lhipp_surface = lhipp_surface,
    lhipp_groups = lhipp_groups,
    lhipp_ids = lhipp_ids,
    lhipp_scans = lhipp_scans,
    lhipp_partition = lhipp_partition
  )
  list2env(out_list, envir = .GlobalEnv)
}
