get_init_params <- function(init_lpme, centers) {
  init_params <- list()
  init_reconstructions <- list()

  for (partition_idx in seq_len(length(init_lpme))) {
    partition_reconstructions <- calculate_lpme_reconstructions(
      init_lpme[[partition_idx]],
      centers[[partition_idx]]
    )

    init_params[[partition_idx]] <- partition_reconstructions$projections
    init_reconstructions[[
      partition_idx
    ]] <- partition_reconstructions$reconstructions
  }

  init_reconstructions <- do.call(rbind, init_reconstructions)

  list(
    params = init_params,
    reconstructions = init_reconstructions
  )
}
