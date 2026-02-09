fit_mnist_additive_model <- function(
  x,
  params,
  weights,
  lambda,
  k,
  ids,
  epsilon = 0.05,
  max_iter = 100
) {
  require(pme, quietly = TRUE)
  require(purrr, quietly = TRUE)

  D <- ncol(x)
  d <- ncol(params)

  fold_vec_full <- rep(1:k, ceiling(length(unique(ids)) / k))
  id_folds <- fold_vec_full[
    sample(
      seq_along(unique(ids)),
      length(unique(ids)),
      replace = FALSE
    )
  ]

  print("Initializing population model")
  init_population_embedding <- fit_weighted_spline(
    x,
    params,
    weights,
    lambda,
    id_folds[ids]
  )

  print("Initializing ID models")
  init_id_embeddings <- list()
  id_x <- list()
  id_params <- list()
  id_weights <- list()
  id_coefs <- list()
  mean_x <- list()

  for (id_idx in seq_along(unique(ids))) {
    id_x[[id_idx]] <- x[ids == unique(ids)[id_idx], ]
    id_params[[id_idx]] <- params[ids == unique(ids)[id_idx], ] %>%
      matrix(ncol = 1)
    # id_weights[[id_idx]] <- weights[ids == unique(ids)[id_idx]]
    id_weights[[id_idx]] <- rep(1 / nrow(id_x[[id_idx]]), nrow(id_x[[id_idx]]))
    fold_vec <- rep(1:k, ceiling(nrow(id_x[[id_idx]]) / k))
    id_folds <- fold_vec[
      sample(
        seq_len(nrow(id_x[[id_idx]])),
        size = nrow(id_x[[id_idx]]),
        replace = FALSE
      )
    ]

    id_preds <- map(
      seq_len(nrow(id_params[[id_idx]])),
      ~ init_population_embedding$embedding_map(id_params[[id_idx]][.x, ])
    ) %>%
      reduce(rbind)
    init_id_embeddings[[id_idx]] <- fit_weighted_spline(
      id_x[[id_idx]] - id_preds,
      id_params[[id_idx]],
      id_weights[[id_idx]],
      lambda,
      id_folds
    )

    mean_x[[id_idx]] <- map(
      seq_len(nrow(id_params[[id_idx]])),
      ~ init_id_embeddings[[id_idx]]$embedding_map(id_params[[id_idx]][.x, ])
    ) %>%
      reduce(rbind) %>%
      colMeans()
  }

  epsilon_hat <- 2 * epsilon
  n <- 0
  population_embedding <- init_population_embedding
  id_embeddings <- init_id_embeddings

  x_preds <- map(
    seq_len(nrow(x)),
    ~ population_embedding$embedding_map(params[.x]) +
      id_embeddings[[ids[.x]]]$embedding_map(params[.x])
  ) %>%
    reduce(rbind)

  mse <- map(
    seq_len(nrow(x)),
    ~ weights[.x] * dist_euclidean(x[.x, ], x_preds[.x, ])^2
  ) %>%
    reduce(c) %>%
    mean()

  while ((epsilon_hat > epsilon) & (n <= max_iter)) {
    print(paste0("Starting iteration ", as.character(n + 1)))
    population_embedding_old <- population_embedding
    id_embeddings_old <- id_embeddings
    coef_diff <- 0
    mse_old <- mse

    x_id_preds <- map(
      seq_len(nrow(x)),
      # ~ id_embeddings_old[[ids[.x]]]$embedding_map(params[.x, ]) - mean_x[[id_idx]]
      ~ id_embeddings_old[[ids[.x]]]$embedding_map(params[.x, ])
    ) %>%
      reduce(rbind)

    fold_vec_full <- rep(1:k, ceiling(length(unique(ids)) / k))
    id_folds <- fold_vec_full[
      sample(
        seq_len(length(unique(ids))),
        length(unique(ids)),
        replace = FALSE
      )
    ]

    population_embedding <- fit_weighted_spline(
      x - x_id_preds,
      params,
      weights,
      lambda,
      id_folds[ids]
    )

    coef_diff <- max(
      coef_diff,
      dist_euclidean(
        population_embedding$coefs,
        population_embedding_old$coefs
      ) /
        norm_euclidean(population_embedding_old$coefs)
    )

    for (id_idx in seq_along(length(unique(ids)))) {
      id_x[[id_idx]] <- x[ids == unique(ids)[id_idx], ]
      id_params[[id_idx]] <- params[ids == unique(ids)[id_idx], ] %>%
        matrix(ncol = 1)
      # id_weights[[id_idx]] <- weights[ids == unique(ids)[id_idx]]
      id_weights[[id_idx]] <- rep(
        1 / nrow(id_x[[id_idx]]),
        nrow(id_x[[id_idx]])
      )

      fold_vec <- rep(1:k, ceiling(nrow(id_x[[id_idx]]) / k))
      id_folds <- fold_vec[sample(
        seq_len(nrow(id_x[[id_idx]])),
        size = nrow(id_x[[id_idx]]),
        replace = FALSE
      )]

      id_preds <- map(
        seq_len(nrow(id_params[[id_idx]])),
        ~ population_embedding$embedding_map(id_params[[id_idx]][.x, ])
      ) %>%
        reduce(rbind)
      id_embeddings[[id_idx]] <- fit_weighted_spline(
        id_x[[id_idx]] - id_preds,
        id_params[[id_idx]],
        id_weights[[id_idx]],
        lambda,
        id_folds
      )

      coef_diff <- max(
        coef_diff,
        dist_euclidean(
          id_embeddings[[id_idx]]$coefs,
          id_embeddings_old[[id_idx]]$coefs
        ) /
          norm_euclidean(id_embeddings_old[[id_idx]]$coefs)
      )

      mean_x[[id_idx]] <- map(
        seq_len(nrow(id_params[[id_idx]])),
        ~ id_embeddings[[id_idx]]$embedding_map(id_params[[id_idx]][.x, ])
      ) %>%
        reduce(rbind) %>%
        colMeans()
    }

    n <- n + 1
    plot_id_pred(
      sample(unique(ids), 1),
      population_embedding,
      id_embeddings,
      params,
      ids
    )

    x_preds <- map(
      seq_len(nrow(x)),
      ~ population_embedding$embedding_map(params[.x]) +
        id_embeddings[[ids[.x]]]$embedding_map(params[.x])
    ) %>%
      reduce(rbind)

    mse <- map(
      seq_len(nrow(x)),
      ~ weights[.x] * dist_euclidean(x[.x, ], x_preds[.x, ])^2
    ) %>%
      reduce(c) %>%
      mean()

    mse_ratio <- abs(mse - mse_old) / mse_old

    epsilon_hat <- mse_ratio
    print(
      paste0(
        "MSE Ratio: ",
        as.character(round(epsilon_hat, 3))
      )
    )
    print(
      paste0(
        "Estimated mean squared error: ",
        as.character(round(mse, 5)),
        "; Relative change in mean squared error: ",
        as.character(round(mse_ratio, 5))
      )
    )
  }

  full_embeddings <- list()
  for (id_idx in seq_along(length(unique(ids)))) {
    full_embeddings[[id_idx]] <- function(parameters) {
      population_embedding$embedding_map(parameters) +
        id_embeddings[[id_idx]]$embedding_map(parameters) -
        mean_x[[id_idx]]
    }
  }

  return(
    list(
      embeddings = full_embeddings,
      population_embedding = population_embedding,
      id_embeddings = id_embeddings
    )
  )
}
