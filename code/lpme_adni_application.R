library(foreach)
library(lubridate)
library(plotly)
library(plot3D)
library(pme)
library(pracma)
library(RColorBrewer)
library(Rfast)
library(tidyverse)

source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est_params.R")
source("code/functions/estimate_volume.R")
source("code/functions/interior_identification.R")
source("code/functions/create_cross_section_matrix.R")



lhipp_surface <- read_csv("data/adni_fsl_lhipp_surface.csv")
lthal_surface <- read_csv("data/adni_fsl_lthal_surface.csv")
rhipp_surface <- read_csv("data/adni_fsl_rhipp_surface.csv")
rthal_surface <- read_csv("data/adni_fsl_rthal_surface.csv")

hipp_info <- read_csv("data/adni_fsl_hipp_info.csv")
thal_info <- read_csv("data/adni_fsl_thal_info.csv")

est_hipp_info <- data.frame(
  patno = character(),
  date = numeric(),
  lhipp_data_vol2 = numeric(),
  lhipp_vol_lpme1 = numeric(),
  lhipp_vol_lpme2 = numeric(),
  lhipp_vol_pme1 = numeric(),
  lhipp_vol_pme2 = numeric(),
  rhipp_data_vol2 = numeric(),
  rhipp_vol_lpme1 = numeric(),
  rhipp_vol_lpme2 = numeric(),
  rhipp_vol_pme1 = numeric(),
  rhipp_vol_pme2 = numeric()
)

est_thal_info <- data.frame(
  patno = character(),
  date = numeric(),
  lthal_data_vol2 = numeric(),
  lthal_vol_lpme1 = numeric(),
  lthal_vol_lpme2 = numeric(),
  lthal_vol_pme1 = numeric(),
  lthal_vol_pme2 = numeric(),
  rthal_data_vol2 = numeric(),
  rthal_vol_lpme1 = numeric(),
  rthal_vol_lpme2 = numeric(),
  rthal_vol_pme1 = numeric(),
  rthal_vol_pme2 = numeric()
)

# process scan dates
lhipp_surface <- lhipp_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )
lthal_surface <- lthal_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )
rhipp_surface <- rhipp_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )
rthal_surface <- rthal_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )

hipp_info <- hipp_info %>%
  mutate(
    date = gsub("\\.0", "", date),
    date = ymd_hms(date),
    date = decimal_date(date)
  )
thal_info <- thal_info %>%
  mutate(
    date = gsub("\\.0", "", date),
    date = ymd_hms(date),
    date = decimal_date(date)
  )

lhipp_surface_centers <- lhipp_surface %>%
  group_by(patno, scan_date) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z)),
    sd_x = sd(x),
    sd_y = sd(y),
    sd_z = sd(z)
  )

lthal_surface_centers <- lthal_surface %>%
  group_by(patno, scan_date) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z)),
    sd_x = sd(x),
    sd_y = sd(y),
    sd_z = sd(z)
  )

rhipp_surface_centers <- rhipp_surface %>%
  group_by(patno, scan_date) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z)),
    sd_x = sd(x),
    sd_y = sd(y),
    sd_z = sd(z)
  )

rthal_surface_centers <- rthal_surface %>%
  group_by(patno, scan_date) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z)),
    sd_x = sd(x),
    sd_y = sd(y),
    sd_z = sd(z)
  )

lhipp_bl <- lhipp_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)

lthal_bl <- lthal_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)

rhipp_bl <- rhipp_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)

rthal_bl <- rthal_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)




# center and standardize data
lhipp_surface <- lhipp_surface %>%
  full_join(lhipp_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(lhipp_bl, by = "patno") %>%
  mutate(
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z),
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    # time_from_bl = (scan_date - time_bl) / duration
    time_from_bl = scan_date - time_bl
  )

lthal_surface <- lthal_surface %>%
  full_join(lthal_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(lthal_bl, by = "patno") %>%
  mutate(
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z),
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    # time_from_bl = (scan_date - time_bl) / duration
    time_from_bl = scan_date - time_bl
  )

rhipp_surface <- rhipp_surface %>%
  full_join(rhipp_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(rhipp_bl, by = "patno") %>%
  mutate(
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z),
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    # time_from_bl = (scan_date - time_bl) / duration
    time_from_bl = scan_date - time_bl
  )

rthal_surface <- rthal_surface %>%
  full_join(rthal_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(rthal_bl, by = "patno") %>%
  mutate(
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z),
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    # time_from_bl = (scan_date - time_bl) / duration
    time_from_bl = scan_date - time_bl
  )

lhipp_patnos <- lhipp_surface %>% 
  filter(duration > 2) %>% 
  .$patno %>% 
  unique()

lhipp_surface <- lhipp_surface %>%
  filter(patno %in% lhipp_patnos)

lthal_surface <- lthal_surface %>%
  filter(patno %in% lhipp_patnos)

rhipp_surface <- rhipp_surface %>%
  filter(patno %in% lhipp_patnos)

rthal_surface <- rthal_surface %>%
  filter(patno %in% lhipp_patnos)


lhipp_surface_spherical <- lhipp_surface %>%
  dplyr::select(x, y, z) %>%
  as.matrix() %>%
  cart2sph() %>%
  as_tibble()

lthal_surface_spherical <- lthal_surface %>%
  dplyr::select(x, y, z) %>%
  as.matrix() %>%
  cart2sph() %>%
  as_tibble()

rhipp_surface_spherical <- rhipp_surface %>%
  dplyr::select(x, y, z) %>%
  as.matrix() %>%
  cart2sph() %>%
  as_tibble()

rthal_surface_spherical <- rthal_surface %>%
  dplyr::select(x, y, z) %>%
  as.matrix() %>%
  cart2sph() %>%
  as_tibble()


lhipp_surface <- bind_cols(lhipp_surface, lhipp_surface_spherical)
lthal_surface <- bind_cols(lthal_surface, lthal_surface_spherical)
rhipp_surface <- bind_cols(rhipp_surface, rhipp_surface_spherical)
rthal_surface <- bind_cols(rthal_surface, rthal_surface_spherical)

patnos <- lhipp_surface$patno %>%
  unique()

ncores <- parallel::detectCores()
cl <- parallel::makeCluster(ncores / 2, type = "FORK")
doParallel::registerDoParallel(cl)

set.seed(10283)
foreach(
  patno_val = patnos,
  .packages = c(
    "pme", 
    "plotly", 
    "plot3D", 
    "pracma", 
    "RColorBrewer", 
    "Rfast", 
    "tidyverse"
  ),
  .export = c(
    "calc_pme_est", 
    "calc_lpme_est_params", 
    "estimate_volume", 
    "interior_identification", 
    "create_cross_section_matrix"
  )
) %do% {
  print(patno_val)
  lhipp <- lhipp_surface %>%
    filter(patno == patno_val)
  lthal <- lthal_surface %>%
    filter(patno == patno_val)
  rhipp <- rhipp_surface %>%
    filter(patno == patno_val)
  rthal <- rthal_surface %>%
    filter(patno == patno_val)

  lhipp_mat <- lhipp %>%
    dplyr::select(
      time_from_bl,
      x,
      y,
      z,
      theta,
      phi,
      r
    ) %>%
    as.matrix()

  lthal_mat <- lthal %>%
    filter(patno == patno_val) %>%
    dplyr::select(
      time_from_bl,
      x,
      y,
      z,
      theta,
      phi,
      r
    ) %>%
    as.matrix()

  rhipp_mat <- rhipp %>%
    filter(patno == patno_val) %>%
    dplyr::select(
      time_from_bl,
      x,
      y,
      z,
      theta,
      phi,
      r
    ) %>%
    as.matrix()
  rthal_mat <- rthal %>%
    filter(patno == patno_val) %>%
    dplyr::select(
      time_from_bl,
      x,
      y,
      z,
      theta,
      phi,
      r
    ) %>%
    as.matrix()


  time_vals <- unique(lhipp_mat[, 1])

  print("Fitting LHIPP Models")

  lhipp_lpme_result_isomap <- lpme(
    lhipp_mat,
    d = 2,
    gamma = exp(-15:5),
    lambda = exp(-12:5),
    initialization_algorithm = "isomap",
    init_type = "subsampling",
    verbose = FALSE,
    print_plots = FALSE
  )
  lhipp_lpme_isomap_est <- calc_lpme_est_params(lhipp_lpme_result_isomap, lhipp_mat)
  lhipp_lpme_isomap_vals <- lhipp_lpme_isomap_est$results
  lhipp_lpme_isomap_params <- lhipp_lpme_isomap_est$params

  lhipp_pme_result <- list()
  lhipp_pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- lhipp_mat[lhipp_mat[, 1] == time_vals[t], -1]
    lhipp_pme_result[[t]] <- pme(temp_data, d = 2, lambda = exp(-12:5), verbose = FALSE)
    lhipp_pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(lhipp_pme_result[[t]], temp_data))
  }
  lhipp_pme_vals <- purrr::reduce(lhipp_pme_vals, rbind)

  patno_path <- paste0("results/", patno_val)
  if (!dir.exists(patno_path)) {
    dir.create(patno_path)
  }
  saveRDS(lhipp_lpme_result_isomap, paste0(patno_path, "/adni_lhipp_lpme_isomap.rds"))
  saveRDS(lhipp_pme_result, paste0(patno_path, "/adni_lhipp_pme.rds"))
  write.csv(lhipp_mat, paste0(patno_path, "/adni_lhipp_mat.csv"))
  write.csv(lhipp_lpme_isomap_vals, paste0(patno_path, "/adni_lhipp_lpme_isomap_vals.csv"))
  write.csv(lhipp_pme_vals, paste0(patno_path, "/adni_lhipp_pme_vals.csv"))

  print("Files saved")

  #### VOLUME ESTIMATION

  # Voxel dimension: 1.2mm x 0.9375mm x 0.9375mm
  # Voxel volume: 1.054688 mm^3

  print("LHIPP Volume Estimation 1")

  
  lhipp_rescaled <- lhipp %>%
    mutate(
      # creating centered coordinates at full scale
      x_rescaled = round((x * max_x)),
      y_rescaled = round((y * max_y)),
      z_rescaled = round((z * max_z))
    )
  candidate_x <- seq(min(lhipp_rescaled$x_rescaled), max(lhipp_rescaled$x_rescaled), 1)
  candidate_y <- seq(min(lhipp_rescaled$y_rescaled), max(lhipp_rescaled$y_rescaled), 1)
  candidate_z <- seq(min(lhipp_rescaled$z_rescaled), max(lhipp_rescaled$z_rescaled), 1)

  candidate_voxels <- expand_grid(candidate_x, candidate_y, candidate_z)

  lhipp_pme_volumes <- vector(mode = "numeric", length = length(time_vals))
  lhipp_lpme_volumes <- vector(mode = "numeric", length = length(time_vals))
  # lhipp_pme_interior_plots <- list()
  # lhipp_lpme_interior_plots <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_max_x <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_candidates <- candidate_voxels %>%
      mutate(
        x_scaled = candidate_x / temp_max_x,
        y_scaled = candidate_y / temp_max_y,
        z_scaled = candidate_z / temp_max_z
      )

    temp_pme <- lhipp_pme_result[[time_idx]]

    temp_pme_embedding <- temp_pme$embedding_map
    temp_pme_coefs <- temp_pme$coefs[[which.min(temp_pme$MSD)]]
    temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

    temp_lpme_coefs <- lhipp_lpme_result_isomap$sol_coef_functions[[which.min(lhipp_lpme_result_isomap$msd)]](time_vals[time_idx])
    temp_lpme_params <- lhipp_lpme_result_isomap$parameterization_list[[which.min(lhipp_lpme_result_isomap$msd)]]
    temp_lpme_params <- temp_lpme_params[temp_lpme_params[, 1] == time_vals[time_idx], -1]
    lpme_n_knots <- nrow(temp_lpme_params)
    d <- ncol(temp_lpme_params)
    coef_mat <- matrix(
      temp_lpme_coefs,
      lpme_n_knots + d + 1,
      byrow = TRUE
    )

    temp_lpme_embedding <- function(r) {
      t(coef_mat[1:lpme_n_knots, ]) %*% pme::etaFunc(r, temp_lpme_params, 4 - d) +
        t(coef_mat[(lpme_n_knots + 1):(lpme_n_knots + d + 1), ]) %*% matrix(c(1, r), ncol = 1)
    }

    temp_candidates_red <- temp_candidates[, 4:6]
    lhipp_interior_voxel_pme <- interior_identification(
      temp_pme_embedding, 
      temp_pme_coefs, 
      temp_pme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )
    lhipp_interior_voxel_lpme <- interior_identification(
      temp_lpme_embedding, 
      coef_mat, 
      temp_lpme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )
   
    lhipp_interior_points_pme <- temp_candidates[lhipp_interior_voxel_pme, ]
    lhipp_interior_points_lpme <- temp_candidates[lhipp_interior_voxel_lpme, ]

    # lhipp_pme_interior_plots[[time_idx]] <- p
    # lhipp_lpme_interior_plots[[time_idx]] <- p_lpme
    lhipp_pme_volumes[time_idx] <- 1.054688 * nrow(lhipp_interior_points_pme)
    lhipp_lpme_volumes[time_idx] <- 1.054688 * nrow(lhipp_interior_points_lpme)
  }

  ### VOLUME ESTIMATION 2

  print("LHIPP Volume Estimation 2")

  lhipp_lpme_vals <- lhipp_lpme_isomap_vals

  lhipp_data_rescaled <- list()
  lhipp_lpme_rescaled <- list()
  lhipp_pme_rescaled <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_data <- lhipp_mat[lhipp_mat[, 1] == time_vals[time_idx], ]
    temp_lpme <- lhipp_lpme_vals[lhipp_lpme_vals[, 1] == time_vals[time_idx], ]
    temp_pme <- lhipp_pme_vals[lhipp_pme_vals[, 1] == time_vals[time_idx], ]

    temp_max_x <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_data_rescaled <- temp_data[, -(5:7)]
    temp_lpme_rescaled <- temp_lpme[, -(5:7)]
    temp_pme_rescaled <- temp_pme[, -(5:7)]

    temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_max_x) / 1.2)
    temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_max_z) / 0.9375)

    lhipp_data_rescaled[[time_idx]] <- temp_data_rescaled
    lhipp_lpme_rescaled[[time_idx]] <- temp_lpme_rescaled
    lhipp_pme_rescaled[[time_idx]] <- temp_pme_rescaled
  }

  lhipp_data_rescaled_full <- reduce(lhipp_data_rescaled, rbind)
  lhipp_lpme_rescaled_full <- reduce(lhipp_lpme_rescaled, rbind)
  lhipp_pme_rescaled_full <- reduce(lhipp_pme_rescaled, rbind)

  voxel_vol <- 1.2 * 0.9375 * 0.9375

  lhipp_data_volumes2 <- map(
    lhipp_data_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  lhipp_data_volume_points2 <- map(
    lhipp_data_volumes2,
    ~ .x[[2]]
  )

  lhipp_data_volumes2 <- map(
    lhipp_data_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  lhipp_lpme_volumes2 <- map(
    lhipp_lpme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  lhipp_lpme_volume_points2 <- map(
    lhipp_lpme_volumes2,
    ~ .x[[2]]
  )

  lhipp_lpme_volumes2 <- map(
    lhipp_lpme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  lhipp_pme_volumes2 <- map(
    lhipp_pme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  lhipp_pme_volume_points2 <- map(
    lhipp_pme_volumes2,
    ~ .x[[2]]
  )

  lhipp_pme_volumes2 <- map(
    lhipp_pme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  lhipp_data_volume_points_full <- map(
    1:length(lhipp_data_volume_points2),
    ~ cbind(time_vals[.x], lhipp_data_volume_points2[[.x]])
  ) %>%
    reduce(rbind)


  lhipp_lpme_volume_points_full <- map(
    1:length(lhipp_lpme_volume_points2),
    ~ cbind(time_vals[.x], lhipp_lpme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  lhipp_pme_volume_points_full <- map(
    1:length(lhipp_pme_volume_points2),
    ~ cbind(time_vals[.x], lhipp_pme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  print("Fitting LTHAL Models")

  lthal_lpme_result_isomap <- lpme(
    lthal_mat,
    d = 2,
    gamma = exp(-15:5),
    lambda = exp(-12:5),
    initialization_algorithm = "isomap",
    init_type = "subsampling",
    verbose = FALSE,
    print_plots = FALSE
  )
  lthal_lpme_isomap_est <- calc_lpme_est_params(lthal_lpme_result_isomap, lthal_mat)
  lthal_lpme_isomap_vals <- lthal_lpme_isomap_est$results
  lthal_lpme_isomap_params <- lthal_lpme_isomap_est$params

  lthal_pme_result <- list()
  lthal_pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- lthal_mat[lthal_mat[, 1] == time_vals[t], -1]
    lthal_pme_result[[t]] <- pme(temp_data, d = 2, lambda = exp(-12:5), verbose = FALSE)
    lthal_pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(lthal_pme_result[[t]], temp_data))
  }
  lthal_pme_vals <- reduce(lthal_pme_vals, rbind)

  saveRDS(lthal_lpme_result_isomap, paste0(patno_path, "/adni_lthal_lpme_isomap.rds"))
  saveRDS(lthal_pme_result, paste0(patno_path, "/adni_lthal_pme.rds"))
  write.csv(lthal_mat, paste0(patno_path, "/adni_lthal_mat.csv"))
  write.csv(lthal_lpme_isomap_vals, paste0(patno_path, "/adni_lthal_lpme_isomap_vals.csv"))
  write.csv(lthal_pme_vals, paste0(patno_path, "/adni_lthal_pme_vals.csv"))

  print("LTHAL Files Saved")


  print("LTHAL Volume Estimation 1")
  #### VOLUME ESTIMATION

  # Voxel dimension: 1.2mm x 0.9375mm x 0.9375mm
  # Voxel volume: 1.054688 mm^3

  lthal_rescaled <- lthal %>%
    mutate(
      x_rescaled = round((x * max_x)),
      y_rescaled = round((y * max_y)),
      z_rescaled = round((z * max_z))
    )
  candidate_x <- seq(min(lthal_rescaled$x_rescaled), max(lthal_rescaled$x_rescaled), 1)
  candidate_y <- seq(min(lthal_rescaled$y_rescaled), max(lthal_rescaled$y_rescaled), 1)
  candidate_z <- seq(min(lthal_rescaled$z_rescaled), max(lthal_rescaled$z_rescaled), 1)

  candidate_voxels <- expand_grid(candidate_x, candidate_y, candidate_z)

  lthal_pme_volumes <- vector(mode = "numeric", length = length(time_vals))
  lthal_lpme_volumes <- vector(mode = "numeric", length = length(time_vals))
  # lthal_pme_interior_plots <- list()
  # lthal_lpme_interior_plots <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_max_x <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_candidates <- candidate_voxels %>%
      mutate(
        x_scaled = candidate_x / temp_max_x,
        y_scaled = candidate_y / temp_max_y,
        z_scaled = candidate_z / temp_max_z
      )

    temp_pme <- lthal_pme_result[[time_idx]]

    temp_pme_embedding <- temp_pme$embedding_map
    temp_pme_coefs <- temp_pme$coefs[[which.min(temp_pme$MSD)]]
    temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

    temp_lpme_coefs <- lthal_lpme_result_isomap$sol_coef_functions[[which.min(lthal_lpme_result_isomap$msd)]](time_vals[time_idx])
    temp_lpme_params <- lthal_lpme_result_isomap$parameterization_list[[which.min(lthal_lpme_result_isomap$msd)]]
    temp_lpme_params <- temp_lpme_params[temp_lpme_params[, 1] == time_vals[time_idx], -1]
    lpme_n_knots <- nrow(temp_lpme_params)
    d <- ncol(temp_lpme_params)
    coef_mat <- matrix(
      temp_lpme_coefs,
      lpme_n_knots + d + 1,
      byrow = TRUE
    )

    temp_lpme_embedding <- function(r) {
      t(coef_mat[1:lpme_n_knots, ]) %*% pme::etaFunc(r, temp_lpme_params, 4 - d) +
        t(coef_mat[(lpme_n_knots + 1):(lpme_n_knots + d + 1), ]) %*% matrix(c(1, r), ncol = 1)
    }
    
    temp_candidates_red <- temp_candidates[, 4:6]
    lthal_interior_voxel_pme <- interior_identification(
      temp_pme_embedding, 
      temp_pme_coefs, 
      temp_pme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )
    lthal_interior_voxel_lpme <- interior_identification(
      temp_lpme_embedding, 
      coef_mat, 
      temp_lpme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )

    lthal_interior_points_pme <- temp_candidates[lthal_interior_voxel_pme, ]
    lthal_interior_points_lpme <- temp_candidates[lthal_interior_voxel_lpme, ]

    # lthal_pme_interior_plots[[time_idx]] <- p
    # lthal_lpme_interior_plots[[time_idx]] <- p_lpme
    lthal_pme_volumes[time_idx] <- 1.054688 * nrow(lthal_interior_points_pme)
    lthal_lpme_volumes[time_idx] <- 1.054688 * nrow(lthal_interior_points_lpme)
  }

  ### VOLUME ESTIMATION 2

  print("LTHAL Volume Estimation 2")

  lthal_lpme_vals <- lthal_lpme_isomap_vals

  lthal_data_rescaled <- list()
  lthal_lpme_rescaled <- list()
  lthal_pme_rescaled <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_data <- lthal_mat[lthal_mat[, 1] == time_vals[time_idx], ]
    temp_lpme <- lthal_lpme_vals[lthal_lpme_vals[, 1] == time_vals[time_idx], ]
    temp_pme <- lthal_pme_vals[lthal_pme_vals[, 1] == time_vals[time_idx], ]

    temp_max_x <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_data_rescaled <- temp_data[, -(5:7)]
    temp_lpme_rescaled <- temp_lpme[, -(5:7)]
    temp_pme_rescaled <- temp_pme[, -(5:7)]

    temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_max_x) / 1.2)
    temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_max_z) / 0.9375)

    lthal_data_rescaled[[time_idx]] <- temp_data_rescaled
    lthal_lpme_rescaled[[time_idx]] <- temp_lpme_rescaled
    lthal_pme_rescaled[[time_idx]] <- temp_pme_rescaled
  }

  lthal_data_rescaled_full <- reduce(lthal_data_rescaled, rbind)
  lthal_lpme_rescaled_full <- reduce(lthal_lpme_rescaled, rbind)
  lthal_pme_rescaled_full <- reduce(lthal_pme_rescaled, rbind)

  voxel_vol <- 1.2 * 0.9375 * 0.9375

  lthal_data_volumes2 <- map(
    lthal_data_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  lthal_data_volume_points2 <- map(
    lthal_data_volumes2,
    ~ .x[[2]]
  )

  lthal_data_volumes2 <- map(
    lthal_data_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  lthal_lpme_volumes2 <- map(
    lthal_lpme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  lthal_lpme_volume_points2 <- map(
    lthal_lpme_volumes2,
    ~ .x[[2]]
  )

  lthal_lpme_volumes2 <- map(
    lthal_lpme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  lthal_pme_volumes2 <- map(
    lthal_pme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  lthal_pme_volume_points2 <- map(
    lthal_pme_volumes2,
    ~ .x[[2]]
  )

  lthal_pme_volumes2 <- map(
    lthal_pme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  lthal_data_volume_points_full <- map(
    1:length(lthal_data_volume_points2),
    ~ cbind(time_vals[.x], lthal_data_volume_points2[[.x]])
  ) %>%
    reduce(rbind)


  lthal_lpme_volume_points_full <- map(
    1:length(lthal_lpme_volume_points2),
    ~ cbind(time_vals[.x], lthal_lpme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  lthal_pme_volume_points_full <- map(
    1:length(lthal_pme_volume_points2),
    ~ cbind(time_vals[.x], lthal_pme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  print("Fitting RHIPP Models")

  rhipp_lpme_result_isomap <- lpme(
    rhipp_mat,
    d = 2,
    gamma = exp(-15:5),
    lambda = exp(-12:5),
    initialization_algorithm = "isomap",
    init_type = "subsampling",
    verbose = FALSE,
    print_plots = FALSE
  )
  rhipp_lpme_isomap_est <- calc_lpme_est_params(rhipp_lpme_result_isomap, rhipp_mat)
  rhipp_lpme_isomap_vals <- rhipp_lpme_isomap_est$results
  rhipp_lpme_isomap_params <- rhipp_lpme_isomap_est$params

  rhipp_pme_result <- list()
  rhipp_pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- rhipp_mat[rhipp_mat[, 1] == time_vals[t], -1]
    rhipp_pme_result[[t]] <- pme(temp_data, d = 2, lambda = exp(-12:5), verbose = FALSE)
    rhipp_pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(rhipp_pme_result[[t]], temp_data))
  }
  rhipp_pme_vals <- reduce(rhipp_pme_vals, rbind)

  saveRDS(rhipp_lpme_result_isomap, paste0(patno_path, "/adni_rhipp_lpme_isomap.rds"))
  saveRDS(rhipp_pme_result, paste0(patno_path, "/adni_rhipp_pme.rds"))
  write.csv(rhipp_mat, paste0(patno_path, "/adni_rhipp_mat.csv"))
  write.csv(rhipp_lpme_isomap_vals, paste0(patno_path, "/adni_rhipp_lpme_isomap_vals.csv"))
  write.csv(rhipp_pme_vals, paste0(patno_path, "/adni_rhipp_pme_vals.csv"))

  print("RHIPP Files Saved")

  #### VOLUME ESTIMATION

  # Voxel dimension: 1.2mm x 0.9375mm x 0.9375mm
  # Voxel volume: 1.054688 mm^3

  print("RHIPP Volume Estimation 1")

  rhipp_rescaled <- rhipp %>%
    mutate(
      x_rescaled = round((x * max_x)),
      y_rescaled = round((y * max_y)),
      z_rescaled = round((z * max_z))
    )
  candidate_x <- seq(min(rhipp_rescaled$x_rescaled), max(rhipp_rescaled$x_rescaled), 1)
  candidate_y <- seq(min(rhipp_rescaled$y_rescaled), max(rhipp_rescaled$y_rescaled), 1)
  candidate_z <- seq(min(rhipp_rescaled$z_rescaled), max(rhipp_rescaled$z_rescaled), 1)

  candidate_voxels <- expand_grid(candidate_x, candidate_y, candidate_z)

  rhipp_pme_volumes <- vector(mode = "numeric", length = length(time_vals))
  rhipp_lpme_volumes <- vector(mode = "numeric", length = length(time_vals))
  # rhipp_pme_interior_plots <- list()
  # rhipp_lpme_interior_plots <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_max_x <- unique(rhipp[rhipp$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(rhipp[rhipp$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(rhipp[rhipp$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_candidates <- candidate_voxels %>%
      mutate(
        x_scaled = candidate_x / temp_max_x,
        y_scaled = candidate_y / temp_max_y,
        z_scaled = candidate_z / temp_max_z
      )

    temp_pme <- rhipp_pme_result[[time_idx]]

    temp_pme_embedding <- temp_pme$embedding_map
    temp_pme_coefs <- temp_pme$coefs[[which.min(temp_pme$MSD)]]
    temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

    temp_lpme_coefs <- rhipp_lpme_result_isomap$sol_coef_functions[[which.min(rhipp_lpme_result_isomap$msd)]](time_vals[time_idx])
    temp_lpme_params <- rhipp_lpme_result_isomap$parameterization_list[[which.min(rhipp_lpme_result_isomap$msd)]]
    temp_lpme_params <- temp_lpme_params[temp_lpme_params[, 1] == time_vals[time_idx], -1]
    lpme_n_knots <- nrow(temp_lpme_params)
    d <- ncol(temp_lpme_params)
    coef_mat <- matrix(
      temp_lpme_coefs,
      lpme_n_knots + d + 1,
      byrow = TRUE
    )

    temp_lpme_embedding <- function(r) {
      t(coef_mat[1:lpme_n_knots, ]) %*% pme::etaFunc(r, temp_lpme_params, 4 - d) +
        t(coef_mat[(lpme_n_knots + 1):(lpme_n_knots + d + 1), ]) %*% matrix(c(1, r), ncol = 1)
    }
    
    temp_candidates_red <- temp_candidates[, 4:6]
    rhipp_interior_voxel_pme <- interior_identification(
      temp_pme_embedding, 
      temp_pme_coefs, 
      temp_pme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )
    rhipp_interior_voxel_lpme <- interior_identification(
      temp_lpme_embedding, 
      coef_mat, 
      temp_lpme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )

    rhipp_interior_points_pme <- temp_candidates[rhipp_interior_voxel_pme, ]
    rhipp_interior_points_lpme <- temp_candidates[rhipp_interior_voxel_lpme, ]

    # rhipp_pme_interior_plots[[time_idx]] <- p
    # rhipp_lpme_interior_plots[[time_idx]] <- p_lpme
    rhipp_pme_volumes[time_idx] <- 1.054688 * nrow(rhipp_interior_points_pme)
    rhipp_lpme_volumes[time_idx] <- 1.054688 * nrow(rhipp_interior_points_lpme)
  }

  print("RHIPP Volume Estimation 2")
  ### VOLUME ESTIMATION 2

  rhipp_lpme_vals <- rhipp_lpme_isomap_vals

  rhipp_data_rescaled <- list()
  rhipp_lpme_rescaled <- list()
  rhipp_pme_rescaled <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_data <- rhipp_mat[rhipp_mat[, 1] == time_vals[time_idx], ]
    temp_lpme <- rhipp_lpme_vals[rhipp_lpme_vals[, 1] == time_vals[time_idx], ]
    temp_pme <- rhipp_pme_vals[rhipp_pme_vals[, 1] == time_vals[time_idx], ]

    temp_max_x <- unique(rhipp[rhipp$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(rhipp[rhipp$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(rhipp[rhipp$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_data_rescaled <- temp_data[, -(5:7)]
    temp_lpme_rescaled <- temp_lpme[, -(5:7)]
    temp_pme_rescaled <- temp_pme[, -(5:7)]

    temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_max_x) / 1.2)
    temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_max_z) / 0.9375)

    rhipp_data_rescaled[[time_idx]] <- temp_data_rescaled
    rhipp_lpme_rescaled[[time_idx]] <- temp_lpme_rescaled
    rhipp_pme_rescaled[[time_idx]] <- temp_pme_rescaled
  }

  rhipp_data_rescaled_full <- reduce(rhipp_data_rescaled, rbind)
  rhipp_lpme_rescaled_full <- reduce(rhipp_lpme_rescaled, rbind)
  rhipp_pme_rescaled_full <- reduce(rhipp_pme_rescaled, rbind)

  voxel_vol <- 1.2 * 0.9375 * 0.9375

  rhipp_data_volumes2 <- map(
    rhipp_data_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  rhipp_data_volume_points2 <- map(
    rhipp_data_volumes2,
    ~ .x[[2]]
  )

  rhipp_data_volumes2 <- map(
    rhipp_data_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  rhipp_lpme_volumes2 <- map(
    rhipp_lpme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  rhipp_lpme_volume_points2 <- map(
    rhipp_lpme_volumes2,
    ~ .x[[2]]
  )

  rhipp_lpme_volumes2 <- map(
    rhipp_lpme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  rhipp_pme_volumes2 <- map(
    rhipp_pme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  rhipp_pme_volume_points2 <- map(
    rhipp_pme_volumes2,
    ~ .x[[2]]
  )

  rhipp_pme_volumes2 <- map(
    rhipp_pme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  rhipp_data_volume_points_full <- map(
    1:length(rhipp_data_volume_points2),
    ~ cbind(time_vals[.x], rhipp_data_volume_points2[[.x]])
  ) %>%
    reduce(rbind)


  rhipp_lpme_volume_points_full <- map(
    1:length(rhipp_lpme_volume_points2),
    ~ cbind(time_vals[.x], rhipp_lpme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  rhipp_pme_volume_points_full <- map(
    1:length(rhipp_pme_volume_points2),
    ~ cbind(time_vals[.x], rhipp_pme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  print("Fitting RTHAL Models")
  rthal_lpme_result_isomap <- lpme(
    rthal_mat,
    d = 2,
    gamma = exp(-15:5),
    lambda = exp(-12:5),
    initialization_algorithm = "isomap",
    init_type = "subsampling",
    verbose = FALSE,
    print_plots = FALSE
  )
  rthal_lpme_isomap_est <- calc_lpme_est_params(rthal_lpme_result_isomap, rthal_mat)
  rthal_lpme_isomap_vals <- rthal_lpme_isomap_est$results
  rthal_lpme_isomap_params <- rthal_lpme_isomap_est$params

  rthal_pme_result <- list()
  rthal_pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- rthal_mat[rthal_mat[, 1] == time_vals[t], -1]
    rthal_pme_result[[t]] <- pme(temp_data, d = 2, lambda = exp(-12:5), verbose = FALSE)
    rthal_pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(rthal_pme_result[[t]], temp_data))
  }
  rthal_pme_vals <- reduce(rthal_pme_vals, rbind)

  saveRDS(rthal_lpme_result_isomap, paste0(patno_path, "/adni_rthal_lpme_isomap.rds"))
  saveRDS(rthal_pme_result, paste0(patno_path, "/adni_rthal_pme.rds"))
  write.csv(rthal_mat, paste0(patno_path, "/adni_rthal_mat.csv"))
  write.csv(rthal_lpme_isomap_vals, paste0(patno_path, "/adni_rthal_lpme_isomap_vals.csv"))
  write.csv(rthal_pme_vals, paste0(patno_path, "/adni_rthal_pme_vals.csv"))

  print("RTHAL Files Saved")

  #### VOLUME ESTIMATION

  # Voxel dimension: 1.2mm x 0.9375mm x 0.9375mm
  # Voxel volume: 1.054688 mm^3

  print("RTHAL Volume Estimation 1")

  rthal_rescaled <- rthal %>%
    mutate(
      x_rescaled = round((x * max_x)),
      y_rescaled = round((y * max_y)),
      z_rescaled = round((z * max_z))
    )
  candidate_x <- seq(min(rthal_rescaled$x_rescaled), max(rthal_rescaled$x_rescaled), 1)
  candidate_y <- seq(min(rthal_rescaled$y_rescaled), max(rthal_rescaled$y_rescaled), 1)
  candidate_z <- seq(min(rthal_rescaled$z_rescaled), max(rthal_rescaled$z_rescaled), 1)

  candidate_voxels <- expand_grid(candidate_x, candidate_y, candidate_z)

  rthal_pme_volumes <- vector(mode = "numeric", length = length(time_vals))
  rthal_lpme_volumes <- vector(mode = "numeric", length = length(time_vals))
  # rthal_pme_interior_plots <- list()
  # rthal_lpme_interior_plots <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_max_x <- unique(rthal[rthal$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(rthal[rthal$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(rthal[rthal$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_candidates <- candidate_voxels %>%
      mutate(
        x_scaled = candidate_x / temp_max_x,
        y_scaled = candidate_y / temp_max_y,
        z_scaled = candidate_z / temp_max_z
      )

    temp_pme <- rthal_pme_result[[time_idx]]

    temp_pme_embedding <- temp_pme$embedding_map
    temp_pme_coefs <- temp_pme$coefs[[which.min(temp_pme$MSD)]]
    temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

    temp_lpme_coefs <- rthal_lpme_result_isomap$sol_coef_functions[[which.min(rthal_lpme_result_isomap$msd)]](time_vals[time_idx])
    temp_lpme_params <- rthal_lpme_result_isomap$parameterization_list[[which.min(rthal_lpme_result_isomap$msd)]]
    temp_lpme_params <- temp_lpme_params[temp_lpme_params[, 1] == time_vals[time_idx], -1]
    lpme_n_knots <- nrow(temp_lpme_params)
    d <- ncol(temp_lpme_params)
    coef_mat <- matrix(
      temp_lpme_coefs,
      lpme_n_knots + d + 1,
      byrow = TRUE
    )

    temp_lpme_embedding <- function(r) {
      t(coef_mat[1:lpme_n_knots, ]) %*% pme::etaFunc(r, temp_lpme_params, 4 - d) +
        t(coef_mat[(lpme_n_knots + 1):(lpme_n_knots + d + 1), ]) %*% matrix(c(1, r), ncol = 1)
    }

    temp_candidates_red <- temp_candidates[, 4:6]
    rthal_interior_voxel_pme <- interior_identification(
      temp_pme_embedding, 
      temp_pme_coefs, 
      temp_pme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )
    rthal_interior_voxel_lpme <- interior_identification(
      temp_lpme_embedding, 
      coef_mat, 
      temp_lpme_params, 
      temp_candidates_red, 
      c(0, 0, 0)
    )
    
    rthal_interior_points_pme <- temp_candidates[rthal_interior_voxel_pme, ]
    rthal_interior_points_lpme <- temp_candidates[rthal_interior_voxel_lpme, ]

    # rthal_pme_interior_plots[[time_idx]] <- p
    # rthal_lpme_interior_plots[[time_idx]] <- p_lpme
    rthal_pme_volumes[time_idx] <- 1.054688 * nrow(rthal_interior_points_pme)
    rthal_lpme_volumes[time_idx] <- 1.054688 * nrow(rthal_interior_points_lpme)
  }

  ### VOLUME ESTIMATION 2

  print("RTHAL Volume Estimation 2")
  rthal_lpme_vals <- rthal_lpme_isomap_vals

  rthal_data_rescaled <- list()
  rthal_lpme_rescaled <- list()
  rthal_pme_rescaled <- list()

  for (time_idx in 1:length(time_vals)) {
    temp_data <- rthal_mat[rthal_mat[, 1] == time_vals[time_idx], ]
    temp_lpme <- rthal_lpme_vals[rthal_lpme_vals[, 1] == time_vals[time_idx], ]
    temp_pme <- rthal_pme_vals[rthal_pme_vals[, 1] == time_vals[time_idx], ]

    temp_max_x <- unique(rthal[rthal$time_from_bl == time_vals[time_idx], ]$max_x)
    temp_max_y <- unique(rthal[rthal$time_from_bl == time_vals[time_idx], ]$max_y)
    temp_max_z <- unique(rthal[rthal$time_from_bl == time_vals[time_idx], ]$max_z)

    temp_data_rescaled <- temp_data[, -(5:7)]
    temp_lpme_rescaled <- temp_lpme[, -(5:7)]
    temp_pme_rescaled <- temp_pme[, -(5:7)]

    temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_max_x) / 1.2)
    temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_max_z) / 0.9375)

    temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_max_x) / 1.2)
    temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_max_y) / 0.9375)
    temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_max_z) / 0.9375)

    rthal_data_rescaled[[time_idx]] <- temp_data_rescaled
    rthal_lpme_rescaled[[time_idx]] <- temp_lpme_rescaled
    rthal_pme_rescaled[[time_idx]] <- temp_pme_rescaled
  }

  rthal_data_rescaled_full <- reduce(rthal_data_rescaled, rbind)
  rthal_lpme_rescaled_full <- reduce(rthal_lpme_rescaled, rbind)
  rthal_pme_rescaled_full <- reduce(rthal_pme_rescaled, rbind)

  voxel_vol <- 1.2 * 0.9375 * 0.9375

  rthal_data_volumes2 <- map(
    rthal_data_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  rthal_data_volume_points2 <- map(
    rthal_data_volumes2,
    ~ .x[[2]]
  )

  rthal_data_volumes2 <- map(
    rthal_data_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  rthal_lpme_volumes2 <- map(
    rthal_lpme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  rthal_lpme_volume_points2 <- map(
    rthal_lpme_volumes2,
    ~ .x[[2]]
  )

  rthal_lpme_volumes2 <- map(
    rthal_lpme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  rthal_pme_volumes2 <- map(
    rthal_pme_rescaled,
    ~ estimate_volume(.x[, -1], voxel_vol)
  )

  rthal_pme_volume_points2 <- map(
    rthal_pme_volumes2,
    ~ .x[[2]]
  )

  rthal_pme_volumes2 <- map(
    rthal_pme_volumes2,
    ~ .x[[1]]
  ) %>%
    reduce(c)

  rthal_data_volume_points_full <- map(
    1:length(rthal_data_volume_points2),
    ~ cbind(time_vals[.x], rthal_data_volume_points2[[.x]])
  ) %>%
    reduce(rbind)


  rthal_lpme_volume_points_full <- map(
    1:length(rthal_lpme_volume_points2),
    ~ cbind(time_vals[.x], rthal_lpme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  rthal_pme_volume_points_full <- map(
    1:length(rthal_pme_volume_points2),
    ~ cbind(time_vals[.x], rthal_pme_volume_points2[[.x]])
  ) %>%
    reduce(rbind)

  dates <- unique(lhipp$time_bl) + time_vals

  print("Creating Temporary Files...")

  temp_hipp_info <- tibble(
    patno = patno_val,
    date = dates,
    lhipp_data_vol2 = lhipp_data_volumes2,
    lhipp_vol_lpme1 = lhipp_lpme_volumes,
    lhipp_vol_lpme2 = lhipp_lpme_volumes2,
    lhipp_vol_pme1 = lhipp_pme_volumes,
    lhipp_vol_pme2 = lhipp_pme_volumes2,
    rhipp_data_vol2 = rhipp_data_volumes2,
    rhipp_vol_lpme1 = rhipp_lpme_volumes,
    rhipp_vol_lpme2 = rhipp_lpme_volumes2,
    rhipp_vol_pme1 = rhipp_pme_volumes,
    rhipp_vol_pme2 = rhipp_pme_volumes2
  )

  temp_thal_info <- tibble(
    patno = patno_val,
    date = dates,
    lthal_data_vol2 = lthal_data_volumes2,
    lthal_vol_lpme1 = lthal_lpme_volumes,
    lthal_vol_lpme2 = lthal_lpme_volumes2,
    lthal_vol_pme1 = lthal_pme_volumes,
    lthal_vol_pme2 = lthal_pme_volumes2,
    rthal_data_vol2 = rthal_data_volumes2,
    rthal_vol_lpme1 = rthal_lpme_volumes,
    rthal_vol_lpme2 = rthal_lpme_volumes2,
    rthal_vol_pme1 = rthal_pme_volumes,
    rthal_vol_pme2 = rthal_pme_volumes2
  )

  est_hipp_info <- bind_rows(est_hipp_info, temp_hipp_info)
  est_thal_info <- bind_rows(est_thal_info, temp_thal_info)

  write.csv(temp_hipp_info, file = paste0(patno_path, "/hipp_info.csv"))
  write.csv(temp_thal_info, file = paste0(patno_path, "/thal_info.csv"))

  print("Temporary Files Saved")

  p <- plot_ly(
    x = lhipp_mat[, 2],
    y = lhipp_mat[, 3],
    z = lhipp_mat[, 4],
    frame = lhipp_mat[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3),
    name = "Data"
  ) %>%
    add_markers(
      x = lhipp_lpme_isomap_vals[, 2],
      y = lhipp_lpme_isomap_vals[, 3],
      z = lhipp_lpme_isomap_vals[, 4],
      frame = lhipp_lpme_isomap_vals[, 1],
      name = "LPME"
    ) %>%
    add_markers(
      x = lhipp_pme_vals[, 2],
      y = lhipp_pme_vals[, 3],
      z = lhipp_pme_vals[, 4],
      frame = lhipp_pme_vals[, 1],
      name = "PME"
    )

  lhipp_lpme_plots <- list()
  lhipp_pme_plots <- list()
  lhipp_data_plots <- list()

  colors <- brewer.pal(3, "Set1")

  png(
    paste0(patno_path, "/adni_lhipp_data_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) > 1e-10) {
      next
    } else {
     temp_data <- lhipp_mat[lhipp_mat[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_data[, 2],
        y = temp_data[, 3],
        z = temp_data[, 4],
        pch = 20,
        col = alpha(colors[1], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      ) 
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_lhipp_lpme_isomap_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) > 1e-10) {
      next
    } else {
      temp_lpme <- lhipp_lpme_isomap_vals[lhipp_lpme_isomap_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_lpme[, 2],
        y = temp_lpme[, 3],
        z = temp_lpme[, 4],
        pch = 20,
        col = alpha(colors[2], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_lhipp_pme_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) > 1e-10) {
      next
    } else {
      temp_pme <- lhipp_pme_vals[lhipp_pme_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_pme[, 2],
        y = temp_pme[, 3],
        z = temp_pme[, 4],
        pch = 20,
        col = alpha(colors[3], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  lhipp_cross_section <- create_cross_section_matrix(
    list(lhipp_mat, lhipp_lpme_isomap_vals, lhipp_pme_vals),
    4,
    c(0.005, 0.005, 0.005),
    nrow = 2,
    lhipp = lhipp
  )
  ggsave(
    paste0(patno_path, "/adni_lhipp_cross_section.png"),
    plot = lhipp_cross_section,
    height = 5000,
    width = 7000,
    units = "px",
    dpi = 1250
  )

  p <- plot_ly(
    x = lthal_mat[, 2],
    y = lthal_mat[, 3],
    z = lthal_mat[, 4],
    frame = lthal_mat[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3),
    name = "Data"
  ) %>%
    add_markers(
      x = lthal_lpme_isomap_vals[, 2],
      y = lthal_lpme_isomap_vals[, 3],
      z = lthal_lpme_isomap_vals[, 4],
      frame = lthal_lpme_isomap_vals[, 1],
      name = "LPME"
    ) %>%
    add_markers(
      x = lthal_pme_vals[, 2],
      y = lthal_pme_vals[, 3],
      z = lthal_pme_vals[, 4],
      frame = lthal_pme_vals[, 1],
      name = "PME"
    )

  lthal_lpme_plots <- list()
  lthal_pme_plots <- list()
  lthal_data_plots <- list()

  colors <- brewer.pal(3, "Set1")

  png(
    paste0(patno_path, "/adni_lthal_data_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) > 1e-10) {
      next
    } else {
      temp_data <- lthal_mat[lthal_mat[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_data[, 2],
        y = temp_data[, 3],
        z = temp_data[, 4],
        pch = 20,
        col = alpha(colors[1], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_lthal_lpme_isomap_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) > 1e-10) {
      next
    } else {
      temp_lpme <- lthal_lpme_isomap_vals[lthal_lpme_isomap_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_lpme[, 2],
        y = temp_lpme[, 3],
        z = temp_lpme[, 4],
        pch = 20,
        col = alpha(colors[2], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_lthal_pme_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) > 1e-10) {
      next
    } else {
      temp_pme <- lthal_pme_vals[lthal_pme_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_pme[, 2],
        y = temp_pme[, 3],
        z = temp_pme[, 4],
        pch = 20,
        col = alpha(colors[3], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  lthal_cross_section <- create_cross_section_matrix(
    list(lthal_mat, lthal_lpme_isomap_vals, lthal_pme_vals),
    4,
    c(0.005, 0.005, 0.005),
    nrow = 2,
    lhipp = lthal
  )
  ggsave(
    paste0(patno_path, "/adni_lthal_cross_section.png"),
    plot = lthal_cross_section,
    height = 5000,
    width = 7000,
    units = "px",
    dpi = 1250
  )

  p <- plot_ly(
    x = rhipp_mat[, 2],
    y = rhipp_mat[, 3],
    z = rhipp_mat[, 4],
    frame = rhipp_mat[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3),
    name = "Data"
  ) %>%
    add_markers(
      x = rhipp_lpme_isomap_vals[, 2],
      y = rhipp_lpme_isomap_vals[, 3],
      z = rhipp_lpme_isomap_vals[, 4],
      frame = rhipp_lpme_isomap_vals[, 1],
      name = "LPME"
    ) %>%
    add_markers(
      x = rhipp_pme_vals[, 2],
      y = rhipp_pme_vals[, 3],
      z = rhipp_pme_vals[, 4],
      frame = rhipp_pme_vals[, 1],
      name = "PME"
    )

  rhipp_lpme_plots <- list()
  rhipp_pme_plots <- list()
  rhipp_data_plots <- list()

  colors <- brewer.pal(3, "Set1")

  png(
    paste0(patno_path, "/adni_rhipp_data_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) != 0) {
      next
    } else {
      temp_data <- rhipp_mat[rhipp_mat[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_data[, 2],
        y = temp_data[, 3],
        z = temp_data[, 4],
        pch = 20,
        col = alpha(colors[1], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(rhipp$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_rhipp_lpme_isomap_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) != 0) {
      next
    } else {
      temp_lpme <- rhipp_lpme_isomap_vals[rhipp_lpme_isomap_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_lpme[, 2],
        y = temp_lpme[, 3],
        z = temp_lpme[, 4],
        pch = 20,
        col = alpha(colors[2], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(rhipp$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_rhipp_pme_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) != 0) {
      next
    } else {
      temp_pme <- rhipp_pme_vals[rhipp_pme_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_pme[, 2],
        y = temp_pme[, 3],
        z = temp_pme[, 4],
        pch = 20,
        col = alpha(colors[3], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(rhipp$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  rhipp_cross_section <- create_cross_section_matrix(
    list(rhipp_mat, rhipp_lpme_isomap_vals, rhipp_pme_vals),
    4,
    c(0.005, 0.005, 0.005),
    nrow = 2,
    lhipp = rhipp
  )
  ggsave(
    paste0(patno_path, "/adni_rhipp_cross_section.png"),
    plot = rhipp_cross_section,
    height = 5000,
    width = 7000,
    units = "px",
    dpi = 1250
  )

  p <- plot_ly(
    x = rthal_mat[, 2],
    y = rthal_mat[, 3],
    z = rthal_mat[, 4],
    frame = rthal_mat[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3),
    name = "Data"
  ) %>%
    add_markers(
      x = rthal_lpme_isomap_vals[, 2],
      y = rthal_lpme_isomap_vals[, 3],
      z = rthal_lpme_isomap_vals[, 4],
      frame = rthal_lpme_isomap_vals[, 1],
      name = "LPME"
    ) %>%
    add_markers(
      x = rthal_pme_vals[, 2],
      y = rthal_pme_vals[, 3],
      z = rthal_pme_vals[, 4],
      frame = rthal_pme_vals[, 1],
      name = "PME"
    )

  rthal_lpme_plots <- list()
  rthal_pme_plots <- list()
  rthal_data_plots <- list()

  colors <- brewer.pal(3, "Set1")

  png(
    paste0(patno_path, "/adni_rthal_data_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) != 0) {
      next
    } else {
      temp_data <- rthal_mat[rthal_mat[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_data[, 2],
        y = temp_data[, 3],
        z = temp_data[, 4],
        pch = 20,
        col = alpha(colors[1], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(rthal$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_rthal_lpme_isomap_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) != 0) {
      next
    } else {
      temp_lpme <- rthal_lpme_isomap_vals[rthal_lpme_isomap_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_lpme[, 2],
        y = temp_lpme[, 3],
        z = temp_lpme[, 4],
        pch = 20,
        col = alpha(colors[2], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(rthal$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  png(
    paste0(patno_path, "/adni_rthal_pme_plot.png"),
    res = 1000,
    height = 15000,
    width = 3500
  )
  par(oma = c(4, 1, 1, 1), mfrow = c(3, 1), mar = c(2, 2, 1, 1))

  for (t in 1:length(time_vals)) {
    if (mod(t, 2) != 0) {
      next
    } else {
      temp_pme <- rthal_pme_vals[rthal_pme_vals[, 1] == time_vals[t], ]
      scatter3D(
        x = temp_pme[, 2],
        y = temp_pme[, 3],
        z = temp_pme[, 4],
        pch = 20,
        col = alpha(colors[3], 0.5),
        theta = 220,
        ticktype = "detailed",
        main = paste0("Time = ", round((time_vals * max(rthal$duration))[t], 2)),
        xlim = c(-0.15, 0.15),
        ylim = c(-0.15, 0.15),
        zlim = c(-0.15, 0.15)
      )
    }
  }
  dev.off()

  rthal_cross_section <- create_cross_section_matrix(
    list(rthal_mat, rthal_lpme_isomap_vals, rthal_pme_vals),
    4,
    c(0.005, 0.005, 0.005),
    nrow = 2,
    lhipp = rthal
  )
  ggsave(
    paste0(patno_path, "/adni_rthal_cross_section.png"),
    plot = rthal_cross_section,
    height = 5000,
    width = 7000,
    units = "px",
    dpi = 1250
  )
}

parallel::stopCluster(cl)

write.csv(est_hipp_info, "results/est_hipp_info.csv")
write.csv(est_thal_info, "results/est_thal_info.csv")
