library(oro.dicom)
library(oro.nifti)
library(neurobase)
library(fslr)
# library(extrantsr)
library(stringr)
library(foreach)
library(tidyverse)
library(lubridate)
# library(ANTsRCore)

sub_dirs <- list.dirs("data/adni", recursive = FALSE)

nii_scans <- list()

processed_scans <- list()

hipp_mask <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz")
)

mni_nifti <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm.nii.gz")
)

ncores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(ncores, type = "FORK")
doParallel::registerDoParallel(cl = cl)

# foreach (dir_idx = 1:length(sub_dirs)) %dopar% {
foreach(dir_idx = 1:8) %dopar% {
  img_dirs <- list.files(sub_dirs[dir_idx], recursive = TRUE, full.names = TRUE) %>%
    gsub(pattern = "([^/]+$)", replacement = "") %>%
    unique()

  nii_scans[[dir_idx]] <- list()

  for (img_idx in 1:length(img_dirs)) {
    proc_dir <- paste0("data/adni_processed_fsl", gsub(pattern = "data/adni", replacement = "", img_dirs[img_idx]))
    if (file.exists(paste0(proc_dir, "/_all_fast_origsegs.nii.gz"))) {
      next
    } else {
      all_slices <- readDICOM(img_dirs[img_idx])
      nii <- dicom2nifti(all_slices)

      run_first_all(
        nii,
        oprefix = proc_dir,
        verbose = FALSE,
        opts = "-d"
      )
    }
  }
}

parallel::stopCluster(cl = cl)

test_dicom <- readDICOM(
  list.files(sub_dirs[1], recursive = TRUE, full.names = TRUE)[1]
)

adni_info <- read_csv("data/adni_info.csv")
adni_info <- adni_info  |> 
  rename(
    image_id = "Image Data ID",
    subid = "Subject",
    date = "Acq Date"
  ) |>
  mutate(
    date = mdy(date)
  ) |>
  select(
    subid,
    Group,
    Sex,
    Age,
    image_id,
    date
  ) |>
  arrange(subid, date)

patnos <- list.dirs("data/adni_processed_fsl", recursive = FALSE, full.names = FALSE)
processed_dirs <- list.dirs("data/adni_processed_fsl", recursive = FALSE)

lhipp <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)
rhipp <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)

lthal <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)
rthal <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)

lhipp_surface <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)
rhipp_surface <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)

lthal_surface <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)
rthal_surface <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  x = double(),
  y = double(),
  z = double(),
  intensity = double()
)

image_info <- tibble(
  subid = character(),
  Group = character(),
  Sex = character(),
  Age = double(),
  image_id = character(),
  date = Date(),
  scan_id = double(),
  lhipp_vol = double(),
  rhipp_vol = double(),
  lthal_vol = double(),
  rthal_vol = double()
)


mprage_names <- c("MP-RAGE", "MPRAGE", "MP-RAGE_REPEAT", "MPRAGE_REPEAT")
for (dir_idx in 1:length(processed_dirs)) {
  print(patnos[dir_idx])
  scan_dates <- list.dirs(
    paste0(processed_dirs[dir_idx], "/", mprage_names),
    full.names = FALSE,
    recursive = FALSE
  )
  scan_dates <- ymd_hms(scan_dates)
  # scan_dates <- floor_date(scan_dates, unit = "day")

  scan_dirs <- list.files(
    paste0(processed_dirs[dir_idx], "/", mprage_names),
    recursive = TRUE,
    full.names = TRUE
  ) %>%
    gsub(pattern = "([^/]+$)", replacement = "") %>%
    unique()
  for (scan_idx in 1:length(scan_dirs)) {
    print(scan_dates[scan_idx])
    scan_info <- adni_info |>
      filter(
        subid == patnos[dir_idx],
        date == floor_date(scan_dates[scan_idx], unit = "day")
      )

    scan_matches <- which(floor_date(scan_dates[scan_idx], unit = "day") == floor_date(scan_dates, unit = "day"))
    scan_id  <- which(scan_idx == scan_matches)
    scan_info <- scan_info[scan_id, ]

    if (file.exists(paste0(scan_dirs[scan_idx], "_all_fast_origsegs.nii.gz"))) {
      temp_segs <- readnii(paste0(scan_dirs[scan_idx], "_all_fast_origsegs.nii.gz"))
      seg_vols <- fslstats(
        paste0(scan_dirs[scan_idx], "_all_fast_origsegs.nii.gz"),
        opts = "-V",
        ts = TRUE
      )
      temp_lhipp <- temp_segs[, , , 6]
      temp_lhipp_vol <- seg_vols[6] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      lhipp_idx <- which(temp_lhipp > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(lhipp_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(lhipp_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- lhipp_idx[lhipp_idx[, dim] == unique_vals[i, 1] & lhipp_idx[, dim2] == unique_vals[i, 2], dim3]
              min_vals <- vector()
              max_vals <- vector()
              for (v in vals) {
                if (!((v-1) %in% vals)) {
                  min_vals <- c(min_vals, v)
                } else if (!((v+1) %in% vals)) {
                  max_vals <- c(max_vals, v)
                }
              }
              surface_vals <- c(min_vals, max_vals)
              for (v in surface_vals) {
                surface[which(lhipp_idx[, dim] == unique_vals[i, 1] & lhipp_idx[, dim2] == unique_vals[i, 2] & lhipp_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      lhipp_temp <- cbind(
        scan_id,
        lhipp_idx %*% diag(pixdim(as.nifti(temp_lhipp))[2:4]),
        temp_lhipp[lhipp_idx]
      )
      df_names <- c("scan_id", "x", "y", "z", "intensity")
      lhipp_temp_df <- as_tibble(lhipp_temp, .name_repair = "minimal")
      names(lhipp_temp_df) <- df_names
      lhipp_surface_temp <- cbind(
        scan_id,
        lhipp_idx[surface, ] %*% diag(pixdim(as.nifti(temp_lhipp))[2:4]),
        temp_lhipp[lhipp_idx[surface, ]]
      )
      lhipp_surface_temp_df <- as_tibble(
        lhipp_surface_temp, 
        .name_repair = "minimal"
      )
      names(lhipp_surface_temp_df) <- df_names

      lhipp_temp_df <- bind_cols(scan_info, lhipp_temp_df)
      lhipp_surface_temp_df <- bind_cols(scan_info, lhipp_surface_temp_df)

      lhipp <- bind_rows(lhipp, lhipp_temp_df)
      lhipp_surface <- bind_rows(lhipp_surface, lhipp_surface_temp_df)

      temp_rhipp <- temp_segs[, , , 13]
      
      temp_rhipp_vol <- seg_vols[13] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      rhipp_idx <- which(temp_rhipp > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(rhipp_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(rhipp_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- rhipp_idx[rhipp_idx[, dim] == unique_vals[i, 1] & rhipp_idx[, dim2] == unique_vals[i, 2], dim3]
              min_vals <- vector()
              max_vals <- vector()
              for (v in vals) {
                if (!((v-1) %in% vals)) {
                  min_vals <- c(min_vals, v)
                } else if (!((v+1) %in% vals)) {
                  max_vals <- c(max_vals, v)
                }
              }
              surface_vals <- c(min_vals, max_vals)
              for (v in surface_vals) {
                surface[which(rhipp_idx[, dim] == unique_vals[i, 1] & rhipp_idx[, dim2] == unique_vals[i, 2] & rhipp_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      rhipp_temp <- cbind(
        scan_id,
        rhipp_idx %*% diag(pixdim(as.nifti(temp_rhipp))[2:4]),
        temp_rhipp[rhipp_idx]
      )
      rhipp_surface_temp <- cbind(
        scan_id,
        rhipp_idx[surface, ] %*% diag(pixdim(as.nifti(temp_rhipp))[2:4]),
        temp_rhipp[rhipp_idx[surface, ]]
      )

      rhipp_temp_df <- as_tibble(rhipp_temp, .name_repair = "minimal")
      names(rhipp_temp_df) <- df_names
      rhipp_surface_temp_df <- as_tibble(
        rhipp_surface_temp, 
        .name_repair = "minimal"
      )
      names(rhipp_surface_temp_df) <- df_names

      rhipp_temp_df <- bind_cols(scan_info, rhipp_temp_df)
      rhipp_surface_temp_df <- bind_cols(scan_info, rhipp_surface_temp_df)

      rhipp <- bind_rows(rhipp, rhipp_temp_df)
      rhipp_surface <- bind_rows(rhipp_surface, rhipp_surface_temp_df)

      temp_lthal <- temp_segs[, , , 1]
      temp_lthal_vol <- seg_vols[1] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      lthal_idx <- which(temp_lthal > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(lthal_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(lthal_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- lthal_idx[lthal_idx[, dim] == unique_vals[i, 1] & lthal_idx[, dim2] == unique_vals[i, 2], dim3]
              min_vals <- vector()
              max_vals <- vector()
              for (v in vals) {
                if (!((v-1) %in% vals)) {
                  min_vals <- c(min_vals, v)
                } else if (!((v+1) %in% vals)) {
                  max_vals <- c(max_vals, v)
                }
              }
              surface_vals <- c(min_vals, max_vals)
              for (v in surface_vals) {
                surface[which(lthal_idx[, dim] == unique_vals[i, 1] & lthal_idx[, dim2] == unique_vals[i, 2] & lthal_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      lthal_temp <- cbind(
        scan_id,
        lthal_idx %*% diag(pixdim(as.nifti(temp_lthal))[2:4]),
        temp_lthal[lthal_idx]
      )
      lthal_surface_temp <- cbind(
        scan_id,
        lthal_idx[surface, ] %*% diag(pixdim(as.nifti(temp_lthal))[2:4]),
        temp_lthal[lthal_idx[surface, ]]
      )

      lthal_temp_df <- as_tibble(lthal_temp, .name_repair = "minimal")
      names(lthal_temp_df) <- df_names
      lthal_surface_temp_df <- as_tibble(
        lthal_surface_temp, 
        .name_repair = "minimal"
      )
      names(lthal_surface_temp_df) <- df_names

      lthal_temp_df <- bind_cols(scan_info, lthal_temp_df)
      lthal_surface_temp_df <- bind_cols(scan_info, lthal_surface_temp_df)

      lthal <- bind_rows(lthal, lthal_temp_df)
      lthal_surface <- bind_rows(lthal_surface, lthal_surface_temp_df)


      temp_rthal <- temp_segs[, , , 9]
      temp_rthal_vol <- seg_vols[9] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      rthal_idx <- which(temp_rthal > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(rthal_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(rthal_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- rthal_idx[rthal_idx[, dim] == unique_vals[i, 1] & rthal_idx[, dim2] == unique_vals[i, 2], dim3]
              min_vals <- vector()
              max_vals <- vector()
              for (v in vals) {
                if (!((v-1) %in% vals)) {
                  min_vals <- c(min_vals, v)
                } else if (!((v+1) %in% vals)) {
                  max_vals <- c(max_vals, v)
                }
              }
              surface_vals <- c(min_vals, max_vals)
              for (v in surface_vals) {
                surface[which(rthal_idx[, dim] == unique_vals[i, 1] & rthal_idx[, dim2] == unique_vals[i, 2] & rthal_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      rthal_temp <- cbind(
        scan_id,
        rthal_idx %*% diag(pixdim(as.nifti(temp_rthal))[2:4]),
        temp_rthal[rthal_idx]
      )
      rthal_surface_temp <- cbind(
        scan_id,
        rthal_idx[surface, ] %*% diag(pixdim(as.nifti(temp_rthal))[2:4]),
        temp_rthal[rthal_idx[surface, ]]
      )

      rthal_temp_df <- as_tibble(rthal_temp, .name_repair = "minimal")
      names(rthal_temp_df) <- df_names
      rthal_surface_temp_df <- as_tibble(
        rthal_surface_temp, 
        .name_repair = "minimal"
      )
      names(rthal_surface_temp_df) <- df_names

      rthal_temp_df <- bind_cols(scan_info, rthal_temp_df)
      rthal_surface_temp_df <- bind_cols(scan_info, rthal_surface_temp_df)

      rthal <- bind_rows(rthal, rthal_temp_df)
      rthal_surface <- bind_rows(rthal_surface, rthal_surface_temp_df)

      scan_info <- scan_info |>
        mutate(
          scan_id = scan_id,
          lhipp_vol = temp_lhipp_vol,
          rhipp_vol = temp_rhipp_vol,
          lthal_vol = temp_lthal_vol,
          rthal_vol = temp_rthal_vol
        )
      image_info <- bind_rows(image_info, scan_info)
    }
  }
}

write_csv(image_info, "data/adni_fsl_info.csv")

write_csv(lhipp, "data/adni_fsl_lhipp.csv")
write_csv(rhipp, "data/adni_fsl_rhipp.csv")

write_csv(lthal, "data/adni_fsl_lthal.csv")
write_csv(rthal, "data/adni_fsl_rthal.csv")

write_csv(lhipp_surface, "data/adni_fsl_lhipp_surface.csv")
write_csv(rhipp_surface, "data/adni_fsl_rhipp_surface.csv")

write_csv(lthal_surface, "data/adni_fsl_lthal_surface.csv")
write_csv(rthal_surface, "data/adni_fsl_rthal_surface.csv")
