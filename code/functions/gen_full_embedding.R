gen_full_embedding <- function(
  population_embedding,
  group_embeddings,
  id_embeddings,
  group_idx,
  id_idx
) {
  require(foreach, quietly = TRUE)

  embedding_map <- function(parameters) {
    population_embedding_out <- population_embedding$embedding_map(parameters)
    group_embedding_out <- group_embeddings[[group_idx]]$embedding_map(
      parameters
    )
    id_embedding_out <- id_embeddings[[id_idx]]$embedding_map(parameters)

    c(
      population_embedding_out[1],
      population_embedding_out[-1] +
        group_embedding_out[-1] +
        id_embedding_out[-1]
    )
  }

  spline_coef_map <- function(timeval) {
    population_embedding$spline_coef_map(timeval) +
      group_embeddings[[group_idx]]$spline_coef_map(timeval) +
      id_embeddings[[id_idx]]$spline_coef_map(timeval)
  }

  list(
    embedding_map = embedding_map,
    spline_coef_map = spline_coef_map,
    group_idx = group_idx,
    id_idx = id_idx
  )
}
