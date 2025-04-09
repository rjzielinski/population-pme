newtemp <- "~/tmp/rcall"
dir.create(newtemp, recursive = TRUE)
Sys.setenv(TMPDIR = tools::file_path_as_absolute(newtemp))
unlink(tempdir(), recursive = TRUE)

library(divest)
library(doFuture)
library(foreach)
library(fslr)
# library(oro.dicom)
# library(oro.nifti)
library(neurobase)
library(progressr)
# library(extrantsr)
library(stringr)
library(tidyverse)
library(lubridate)
# library(ANTsRCore)
handlers(global = TRUE)

skip_imgs <- c("I436253", "I59843")

img_dirs <- list.files("data/ADNI", recursive = TRUE, full.names = TRUE) %>%
  gsub(pattern = "([^/]+$)", replacement = "") %>%
  unique()

hipp_mask <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz")
)

mni_nifti <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm.nii.gz")
)

ncores <- parallel::detectCores()
plan(cluster, workers = ncores / 2)

segment_imgs <- function(img_dirs) {
  p <- progressor(along = img_dirs)
  error_list <- foreach (
      img_idx = seq_along(img_dirs), 
      # .packages = c("oro.dicom", "oro.nifti", "neurobase", "fslr", "magrittr"),
      # .export = c("img_dirs"),
      .errorhandling = "pass"
    ) %dofuture% {
# foreach(dir_idx = 1:64) %dopar% {
      img_val <- str_split(img_dirs[img_idx], "/")[[1]][6]
      proc_dir <- paste0(
        "data/adni_processed_fsl", 
        gsub(pattern = "data/ADNI", replacement = "", img_dirs[img_idx])
      )
      if (file.exists(paste0(proc_dir, "/_all_fast_origsegs.nii.gz"))) {
        seg_error <- FALSE
      } else if (img_val %in% skip_imgs) {
        seg_error <- TRUE
      } else {
        # all_slices <- readDICOM(img_dirs[img_idx])
        # nii <- dicom2nifti(all_slices)
        nii <- readDicom(img_dirs[img_idx], interactive = FALSE, verbosity = -1)

        first_out <- run_first_all(
          nii,
          oprefix = proc_dir,
          verbose = FALSE,
          opts = "-d"
        )
        p(sprintf("img: %s", img_val))
        seg_error <- is.null(first_out$segmentation)
      }
      seg_error
    }
  return(error_list)
}

error_list <- segment_imgs(img_dirs)

error_list <- reduce(error_list, c)
segmentation_errors <- img_dirs[error_list]

saveRDS(segmentation_errors, "data/segmentation_errors.rds")
