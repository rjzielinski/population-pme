plot_id_pred <- function(id, pop_embedding, id_embeddings, params, ids) {
  # mnist_data <- mnist_x[mnist_id == id, ]
  mnist_data <- input[id_vals == id, ]
  id_params <- params[ids == id, ]

  means <- list()
  for (id_idx in 1:length(id_embeddings)) {
    means[[id_idx]] <- map(
      1:nrow(matrix(params[ids == id, ])),
      ~ id_embeddings[[id_idx]]$embedding_map(matrix(params[ids == id, ])[.x, ])
    ) %>%
      reduce(rbind) %>%
      colMeans()
  }

  sample_pop_params <- seq(min(params), max(params), length.out = 1000)
  sample_id_params <- seq(min(id_params), max(id_params), length.out = 1000)
  pop_out <- map(1:length(sample_pop_params), ~ pop_embedding$embedding_map(sample_pop_params[.x])) %>%
    reduce(rbind)
  id_out <- map(
    1:length(sample_id_params),
    ~ id_embeddings[[id]]$embedding_map(sample_id_params[.x]) - means[[id]]
  ) %>%
    reduce(rbind)
  plt <- ggplot() +
    geom_jitter(aes(x = mnist_data[, 1], y = mnist_data[, 2])) +
    geom_point(aes(x = pop_out[, 1], y = pop_out[, 2]), color = "red") +
    geom_point(aes(x = pop_out[, 1] + id_out[, 1], y = pop_out[, 2] + id_out[, 2]), color = "blue") +
    xlab("X") +
    ylab("Y")
  print(plt)
}
