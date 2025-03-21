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

sub_dirs <- list.dirs("data/ADNI", recursive = FALSE)

nii_scans <- list()

processed_scans <- list()

hipp_mask <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz")
)

mni_nifti <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm.nii.gz")
)

ncores <- parallel::detectCores() / 2
cl <- parallel::makeCluster(ncores, type = "FORK")
doParallel::registerDoParallel(cl = cl)

foreach (dir_idx = 1:length(sub_dirs), .errorhandling = "pass") %dopar% {
# foreach(dir_idx = 1:64) %dopar% {
  img_dirs <- list.files(sub_dirs[dir_idx], recursive = TRUE, full.names = TRUE) %>%
    gsub(pattern = "([^/]+$)", replacement = "") %>%
    unique()

  nii_scans[[dir_idx]] <- list()

  for (img_idx in 1:length(img_dirs)) {
    proc_dir <- paste0("data/adni_processed_fsl", gsub(pattern = "data/ADNI", replacement = "", img_dirs[img_idx]))
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
