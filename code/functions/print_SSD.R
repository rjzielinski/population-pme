print_SSD <- function(SSD_new, SSD_ratio, count) {
  print(
    paste0(
      "For iteration # ",
      as.character(count),
      ", SSD = ",
      as.character(round(SSD_new, 4)),
      " and SSD_ratio = ",
      as.character(round(SSD_ratio, 4)),
      "."
    )
  )
}
