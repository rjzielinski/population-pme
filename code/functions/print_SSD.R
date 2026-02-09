print_ssd <- function(ssd_new, ssd_ratio, count) {
  print(
    paste0(
      "For iteration # ",
      as.character(count),
      ", ssd = ",
      as.character(round(ssd_new, 4)),
      " and ssd_ratio = ",
      as.character(round(ssd_ratio, 4)),
      "."
    )
  )
}
