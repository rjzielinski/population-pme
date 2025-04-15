# newtemp <- "~/tmp/rcall"
# dir.create(newtemp, recursive = TRUE)
# Sys.setenv(TMPDIR = tools::file_path_as_absolute(newtemp))
# unlink(tempdir(), recursive = TRUE)

library(divest)
library(doFuture)
library(dplyr)
library(foreach)
library(fslr)
library(oro.dicom)
library(oro.nifti)
library(neurobase)
library(progressr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyverse)
library(lubridate)
handlers(global = TRUE)

skip_imgs <- c("I436253", "I59843", "I272379", "I135611", "I48503", "I77171")

img_dirs <- list.files("data/ADNI", recursive = TRUE, full.names = TRUE) %>%
  gsub(pattern = "([^/]+$)", replacement = "") %>%
  unique()

hipp_mask <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz")
)

mni_nifti <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm.nii.gz")
)

adni_info <- read_csv("data/adni_info.csv")
adni_info <- adni_info |>
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


# plan(cluster, workers = ncores)
plan(multisession)

segment_imgs <- function(img_dirs) {
  p <- progressor(along = img_dirs)

  image_info <- foreach(
    img_idx = seq_along(img_dirs),
    .errorhandling = "pass"
  ) %dofuture% {
    img_val <- str_split(img_dirs[img_idx], "/")[[1]][6]
    proc_dir <- paste0(
      "data/adni_processed_fsl",
      gsub(pattern = "data/ADNI", replacement = "", img_dirs[img_idx])
    )
    if (file.exists(paste0(proc_dir, "_all_fast_origsegs.nii.gz"))) {
      NULL
    } else if (img_val %in% skip_imgs) {
      NULL
    } else {
      nii <- readDicom(img_dirs[img_idx], interactive = FALSE, verbosity = -1)

      run_first_all(
        nii,
        oprefix = proc_dir,
        verbose = FALSE,
        opts = "-d"
      )
      gc()
    }

    if (file.exists(paste0(proc_dir, "_all_fast_origsegs.nii.gz"))) {
      temp_segs <- readnii(paste0(proc_dir, "_all_fast_origsegs.nii.gz"))
      seg_vols <- fslstats(
        paste0(proc_dir, "_all_fast_origsegs.nii.gz"),
        opts = "-V",
        ts = TRUE
      )
      lhipp_vol <- seg_vols[6] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()

      rhipp_vol <- seg_vols[13] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()

      lthal_vol <- seg_vols[1] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()

      rthal_vol <- seg_vols[9] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()


      scan_info <- tibble(
        image_id = img_val,
        lhipp_vol = lhipp_vol,
        rhipp_vol = rhipp_vol,
        lthal_vol = lthal_vol,
        rthal_vol = rthal_vol
      )

      if (file.exists(paste0(proc_dir, "lhipp.csv"))) {
        NULL
      } else if (img_val %in% skip_imgs) {
        NULL
      } else {
        lhipp <- temp_segs[, , , 6]
        lhipp_idx <- which(lhipp > 0, arr.ind = TRUE)
        surface <- rep(FALSE, nrow(lhipp_idx))
        for (dim in 1:3) {
          for (dim2 in 1:3) {
            if (dim != dim2) {
              dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
              unique_vals <- unique(lhipp_idx[, c(dim, dim2)])
              for (i in seq_len(nrow(unique_vals))) {
                vals <- lhipp_idx[lhipp_idx[, dim] == unique_vals[i, 1] & lhipp_idx[, dim2] == unique_vals[i, 2], dim3]
                min_vals <- vector()
                max_vals <- vector()
                for (v in vals) {
                  if (!((v - 1) %in% vals)) {
                    min_vals <- c(min_vals, v)
                  } else if (!((v + 1) %in% vals)) {
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
        lhipp_mat <- cbind(
          img_val,
          lhipp_idx %*% diag(pixdim(as.nifti(lhipp))[2:4]),
          lhipp[lhipp_idx]
        )
        df_names <- c("scan_id", "x", "y", "z", "intensity")
        lhipp_df <- as_tibble(lhipp_mat, .name_repair = "minimal")
        names(lhipp_df) <- df_names
        lhipp_surface <- cbind(
          img_val,
          lhipp_idx[surface, ] %*% diag(pixdim(as.nifti(lhipp))[2:4]),
          lhipp[lhipp_idx[surface, ]]
        )
        lhipp_surface_df <- as_tibble(
          lhipp_surface,
          .name_repair = "minimal"
        )
        names(lhipp_surface_df) <- df_names

        write_csv(lhipp_df, paste0(proc_dir, "lhipp.csv"))
        write_csv(lhipp_surface_df, paste0(proc_dir, "lhipp_surface.csv"))
      }

      if (file.exists(paste0(proc_dir, "rhipp.csv"))) {
        NULL
      } else if (img_val %in% skip_imgs) {
        NULL
      } else {
        rhipp <- temp_segs[, , , 13]
        rhipp_idx <- which(rhipp > 0, arr.ind = TRUE)
        surface <- rep(FALSE, nrow(rhipp_idx))
        for (dim in 1:3) {
          for (dim2 in 1:3) {
            if (dim != dim2) {
              dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
              unique_vals <- unique(rhipp_idx[, c(dim, dim2)])
              for (i in seq_len(nrow(unique_vals))) {
                vals <- rhipp_idx[rhipp_idx[, dim] == unique_vals[i, 1] & rhipp_idx[, dim2] == unique_vals[i, 2], dim3]
                min_vals <- vector()
                max_vals <- vector()
                for (v in vals) {
                  if (!((v - 1) %in% vals)) {
                    min_vals <- c(min_vals, v)
                  } else if (!((v + 1) %in% vals)) {
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
        rhipp_mat <- cbind(
          img_val,
          rhipp_idx %*% diag(pixdim(as.nifti(rhipp))[2:4]),
          rhipp[rhipp_idx]
        )
        df_names <- c("scan_id", "x", "y", "z", "intensity")
        rhipp_df <- as_tibble(rhipp_mat, .name_repair = "minimal")
        names(rhipp_df) <- df_names
        rhipp_surface <- cbind(
          img_val,
          rhipp_idx[surface, ] %*% diag(pixdim(as.nifti(rhipp))[2:4]),
          rhipp[rhipp_idx[surface, ]]
        )
        rhipp_surface_df <- as_tibble(
          rhipp_surface,
          .name_repair = "minimal"
        )
        names(rhipp_surface_df) <- df_names

        write_csv(rhipp_df, paste0(proc_dir, "rhipp.csv"))
        write_csv(rhipp_surface_df, paste0(proc_dir, "rhipp_surface.csv"))
      }

      if (file.exists(paste0(proc_dir, "lthal.csv"))) {
        NULL
      } else if (img_val %in% skip_imgs) {
        NULL
      } else {
        lthal <- temp_segs[, , , 1]
        lthal_idx <- which(lthal > 0, arr.ind = TRUE)
        surface <- rep(FALSE, nrow(lthal_idx))
        for (dim in 1:3) {
          for (dim2 in 1:3) {
            if (dim != dim2) {
              dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
              unique_vals <- unique(lthal_idx[, c(dim, dim2)])
              for (i in seq_len(nrow(unique_vals))) {
                vals <- lthal_idx[lthal_idx[, dim] == unique_vals[i, 1] & lthal_idx[, dim2] == unique_vals[i, 2], dim3]
                min_vals <- vector()
                max_vals <- vector()
                for (v in vals) {
                  if (!((v - 1) %in% vals)) {
                    min_vals <- c(min_vals, v)
                  } else if (!((v + 1) %in% vals)) {
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
        lthal_mat <- cbind(
          img_val,
          lthal_idx %*% diag(pixdim(as.nifti(lthal))[2:4]),
          lthal[lthal_idx]
        )
        df_names <- c("scan_id", "x", "y", "z", "intensity")
        lthal_df <- as_tibble(lthal_mat, .name_repair = "minimal")
        names(lthal_df) <- df_names
        lthal_surface <- cbind(
          img_val,
          lthal_idx[surface, ] %*% diag(pixdim(as.nifti(lthal))[2:4]),
          lthal[lthal_idx[surface, ]]
        )
        lthal_surface_df <- as_tibble(
          lthal_surface,
          .name_repair = "minimal"
        )
        names(lthal_surface_df) <- df_names

        write_csv(lthal_df, paste0(proc_dir, "lthal.csv"))
        write_csv(lthal_surface_df, paste0(proc_dir, "lthal_surface.csv"))
      }

      if (file.exists(paste0(proc_dir, "rthal.csv"))) {
        NULL
      } else if (img_val %in% skip_imgs) {
        NULL
      } else {
        rthal <- temp_segs[, , , 9]
        rthal_idx <- which(rthal > 0, arr.ind = TRUE)
        surface <- rep(FALSE, nrow(rthal_idx))
        for (dim in 1:3) {
          for (dim2 in 1:3) {
            if (dim != dim2) {
              dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
              unique_vals <- unique(rthal_idx[, c(dim, dim2)])
              for (i in seq_len(nrow(unique_vals))) {
                vals <- rthal_idx[rthal_idx[, dim] == unique_vals[i, 1] & rthal_idx[, dim2] == unique_vals[i, 2], dim3]
                min_vals <- vector()
                max_vals <- vector()
                for (v in vals) {
                  if (!((v - 1) %in% vals)) {
                    min_vals <- c(min_vals, v)
                  } else if (!((v + 1) %in% vals)) {
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
        rthal_mat <- cbind(
          img_val,
          rthal_idx %*% diag(pixdim(as.nifti(rthal))[2:4]),
          rthal[rthal_idx]
        )
        df_names <- c("scan_id", "x", "y", "z", "intensity")
        rthal_df <- as_tibble(rthal_mat, .name_repair = "minimal")
        names(rthal_df) <- df_names
        rthal_surface <- cbind(
          img_val,
          rthal_idx[surface, ] %*% diag(pixdim(as.nifti(rthal))[2:4]),
          rthal[rthal_idx[surface, ]]
        )
        rthal_surface_df <- as_tibble(
          rthal_surface,
          .name_repair = "minimal"
        )
        names(rthal_surface_df) <- df_names

        write_csv(rthal_df, paste0(proc_dir, "rthal.csv"))
        write_csv(rthal_surface_df, paste0(proc_dir, "rthal_surface.csv"))
      }
    }

    p(sprintf("img: %s", img_val))
    scan_info
  }
}

image_info <- segment_imgs(img_dirs)

info_present <- map(
  image_info,
  ~ !(is.null(.x) | "simpleError" %in% class(.x))
) |>
  reduce(c)
image_info <- image_info[info_present] |>
  reduce(bind_rows)


adni_info <- full_join(
  adni_info,
  image_info,
  by = "image_id"
)

write_csv(adni_info, "data/adni_info_full.csv")

processed_dirs <- gsub(
  pattern = "data/ADNI",
  replacement = "data/adni_processed_fsl",
  img_dirs
)
error_list <- map(
  processed_dirs,
  ~ ifelse(
    file.exists(paste0(.x, "/_all_fast_origsegs.nii.gz")),
    FALSE,
    TRUE
  )
) |>
  reduce(c)
segmentation_errors <- img_dirs[error_list]

saveRDS(segmentation_errors, "data/segmentation_errors.rds")
