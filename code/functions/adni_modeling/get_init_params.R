get_init_params <- function(init_lpme, centers) {
  init_params <- list()

  for (partition_idx in seq_len(length(init_lpme))) {
    init_params[[partition_idx]] <- calculate_lpme_reconstructions(
      init_lpme[[partition_idx]],
      centers[[partition_idx]]
    )$projections
  }

  init_params
}
