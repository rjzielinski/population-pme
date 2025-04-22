library(dplyr)
library(progressr)
library(purrr)
library(readr)

handlers(global = TRUE)

combine_csvs <- function(files) {
  df_list <- map(files, read_csv, show_col_types = FALSE)
  combined_df <- reduce(df_list, bind_rows)
  combined_df
}

file_names <- c(
  "lhipp",
  "lhipp_surface",
  "rhipp",
  "rhipp_surface",
  "lthal",
  "lthal_surface",
  "rthal",
  "rthal_surface"
)

for (file_name in file_names) {
  print(paste("Processing", file_name, "files..."))
  file_vector <- list.files(
    "data/adni_processed_fsl",
    pattern = paste0(file_name, ".csv"),
    full.names = TRUE,
    recursive = TRUE
  )

  full_name <- paste0("data/", file_name, "_fsl.csv")
  command_string1 <- paste(
    "head -n 1", 
    file_vector[1], 
    ">", 
    full_name,
    sep = " "
  ) 

  system(command_string1)

  for (f_name in file_vector) {
    command_string_f <- paste(
      "tail -n+2 -q", 
      f_name, 
      ">>", 
      full_name, 
      sep = " "
    )
    system(command_string_f)
  }
}
