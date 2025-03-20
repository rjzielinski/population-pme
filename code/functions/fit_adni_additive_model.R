fit_adni_additive_model <- function(x, params, weights, lambda, k, groups, ids, scans, times, epsilon = 0.05, max_iter = 100) {
  D <- ncol(x)
  d <- ncol(params)

  groups <- as.factor(groups)
  ids = as.factor(ids)
  scans <- as.factor(scans)


  fold_vec_full <- rep(1:k, ceiling(length(unique(scans)) / k))
  scan_folds <- fold_vec_full[
    sample(
      1:length(unique(scans)),
      length(unique(scans)),
      replace = FALSE
    )
  ]

  print("Initializing population model")
  init_population_embedding <- fit_weighted_spline(
    x,
    params,
    weights,
    lambda,
    scan_folds[scans]
  )

  print("Initializing Group models")
  init_group_embeddings <- list()
  group_x <- list()
  group_params <- list()
  group_weights <- list()
  group_coefs <- list()

  for (group_idx in seq_along(unique(groups))) {
    group_set <- groups == unique(groups)[group_idx]
    group_x[[group_idx]] <- x[group_set,]
    group_params[[group_idx]] <- params[group_set,] 
    group_weights[[group_idx]] <- weights[group_set]
    group_ids <- droplevels(ids[group_set])
    if (length(unique(group_ids)) < k) {
      temp_k <- length(unique(group_ids))
    } else {
      temp_k <- k
    }
    id_fold_vec <- rep(1:temp_k, ceiling(length(unique(group_ids)) / temp_k))
    id_folds <- id_fold_vec[
      sample(
        unique(group_ids),
        length(unique(group_ids)),
        replace = FALSE
      )
    ]
    group_preds <- map(
      1:nrow(group_params[[group_idx]]),
      ~ init_population_embedding$embedding_map(group_params[[group_idx]][.x, ])
    ) %>%
      reduce(rbind)
    init_group_embeddings[[group_idx]] <- fit_weighted_spline(
      group_x[[group_idx]] - group_preds,
      group_params[[group_idx]],
      group_weights[[group_idx]],
      lambda,
      id_folds[group_ids]
    )
  }

  print("Initializing ID models")
  init_id_embeddings <- list()
  id_x <- list()
  id_params <- list()
  id_weights <- list()
  id_coefs <- list()
  for (id_idx in 1:length(unique(ids))) {
    id_set <- ids == unique(ids)[id_idx]
    id_group  <- which(unique(groups) == unique(groups[id_set]))
    id_x[[id_idx]] <- x[id_set, ]
    id_params[[id_idx]] <- params[id_set, ]
    id_weights[[id_idx]] <- weights[id_set]
    id_scans <- droplevels(scans[id_set])
    if (length(unique(id_scans)) < k) {
      temp_k <- length(unique(id_scans))
    } else {
      temp_k <- k
    }
    scan_fold_vec <- rep(1:temp_k, ceiling(length(unique(id_scans)) / temp_k))
    scan_folds <- scan_fold_vec[
      sample(
        unique(id_scans),
        length(unique(id_scans)),
        replace = FALSE
      )
    ]
    scan_fold_vals <- scan_folds[id_scans]
    if (length(unique(id_scans)) == 1) {
      scan_fold_vec <- rep(1:k, ceiling(nrow(id_x[[id_idx]]) / k))
      scan_fold_vals <- scan_fold_vec[
        sample(
          1:nrow(id_x[[id_idx]]),
          nrow(id_x[[id_idx]]),
          replace = FALSE
        )
      ]
    }
    id_preds <- map(
      1:nrow(id_params[[id_idx]]),
      ~ (
       init_population_embedding$embedding_map(id_params[[id_idx]][.x, ]) +
       init_group_embeddings[[id_group]]$embedding_map(id_params[[id_idx]][.x, ])
      )
    ) %>%
      reduce(rbind)
    init_id_embeddings[[id_idx]] <- fit_weighted_spline(
      id_x[[id_idx]] - id_preds,
      id_params[[id_idx]],
      id_weights[[id_idx]],
      lambda,
      scan_fold_vals
    )
  }

  print("Initializing Image models")
  init_img_embeddings <- list()
  img_x <- list()
  img_params <- list()
  img_weights <- list()
  img_coefs <- list()
  mean_x <- list()
  for (img_idx in seq_along(unique(scans))) {
    img_set <- scans == unique(scans)[img_idx]
    img_group  <- which(unique(groups) == unique(groups[img_set]))
    img_id <- which(unique(ids) == unique(ids[img_set]))
    img_x[[img_idx]] <- x[img_set, ]
    img_params[[img_idx]] <- params[img_set, ]
    img_weights[[img_idx]] <- weights[img_set]
    obs_fold_vec <- rep(1:k, ceiling(sum(img_set) / k))
    obs_folds <- obs_fold_vec[
      sample(
        1:sum(img_set),
        sum(img_set),
        replace = FALSE
      )
    ]
    img_preds <- map(
      1:nrow(img_params[[img_idx]]),
      ~ (
       init_population_embedding$embedding_map(img_params[[img_idx]][.x, ]) +
       init_group_embeddings[[img_group]]$embedding_map(img_params[[img_idx]][.x, ]) +
       init_id_embeddings[[img_id]]$embedding_map(img_params[[img_idx]][.x, ])
      )
    ) %>%
      reduce(rbind)
    init_img_embeddings[[img_idx]] <- fit_weighted_spline(
      img_x[[img_idx]] - img_preds,
      img_params[[img_idx]],
      img_weights[[img_idx]],
      lambda,
      obs_folds
    )
    mean_x[[img_idx]] <- map(
      1:nrow(img_params[[img_idx]]),
      ~ (
        init_population_embedding$embedding_map(img_params[[img_idx]][.x, ]) +
        init_group_embeddings[[img_group]]$embedding_map(img_params[[img_idx]][.x, ]) +
        init_id_embeddings[[img_id]]$embedding_map(img_params[[img_idx]][.x, ]) +
        init_img_embeddings[[img_idx]]$embedding_map(img_params[[img_idx]][.x, ])
      )
    ) %>%
      reduce(rbind) %>%
      colMeans()
  }

  epsilon_hat <- 2 * epsilon
  n <- 0
  population_embedding <- init_population_embedding
  group_embeddings <- init_group_embeddings
  id_embeddings <- init_id_embeddings
  img_embeddings <- init_img_embeddings

  x_preds <- map(
    1:nrow(x),
    ~ population_embedding$embedding_map(params[.x, ]) +
      group_embeddings[[groups[.x]]]$embedding_map(params[.x, ]) +
      id_embeddings[[ids[.x]]]$embedding_map(params[.x, ]) +
      img_embeddings[[scans[.x]]]$embedding_map(params[.x, ])
  ) %>%
    reduce(rbind)

  mse <- map(1:nrow(x), ~ weights[.x] * dist_euclidean(x[.x, ], x_preds[.x, ])^2) %>%
    reduce(c) %>%
    mean()

  while ((epsilon_hat > epsilon) & (n <= max_iter)) {
    print(paste0("Starting iteration ", as.character(n + 1)))
    population_embedding_old <- population_embedding
    group_embeddings_old <- group_embeddings
    id_embeddings_old <- id_embeddings
    img_embeddings_old <- img_embeddings
    mse_old <- mse

    x_pred_nopop <- map(
      1:nrow(x),
      ~ group_embeddings[[groups[.x]]]$embedding_map(params[.x, ]) +
        id_embeddings[[ids[.x]]]$embedding_map(params[.x, ]) +
        img_embeddings[[scans[.x]]]$embedding_map(params[.x, ])
    ) %>%
      reduce(rbind)

    fold_vec_full <- rep(1:k, ceiling(length(unique(scans)) / k))
    scan_folds <- fold_vec_full[
      sample(
        1:length(unique(scans)),
        length(unique(scans)),
        replace = FALSE
      )
    ]

    population_embedding <- fit_weighted_spline(
      x - x_pred_nopop,
      params,
      weights,
      lambda,
      scan_folds[scans]
    )

    for (group_idx in seq_along(unique(groups))) {
      group_set <- groups == unique(groups)[group_idx]
      group_x[[group_idx]] <- x[group_set, ]
      group_params[[group_idx]] <- params[group_set, ]
      group_weights[[group_idx]] <- weights[group_set]
      group_ids <- droplevels(ids[group_set])
      if (length(unique(group_ids)) < k) {
        temp_k <- length(unique(group_ids))
      } else {
        temp_k <- k
      }
      id_fold_vec <- rep(1:temp_k, ceiling(length(unique(group_ids)) / temp_k))
      id_folds <- id_fold_vec[
        sample(
          unique(group_ids),
          length(unique(group_ids)),
          replace = FALSE
        )
      ]
      group_preds <- map(
        1:nrow(group_params[[group_idx]]),
        ~ (
          population_embedding$embedding_map(group_params[[group_idx]][.x, ]) +
          id_embeddings[[ids[group_set][.x]]]$embedding_map(group_params[[group_idx]][.x, ]) +
          img_embeddings[[scans[group_set][.x]]]$embedding_map(group_params[[group_idx]][.x, ])
        )
      ) %>%
        reduce(rbind)
      group_embeddings[[group_idx]] <- fit_weighted_spline(
        group_x[[group_idx]] - group_preds,
        group_params[[group_idx]],
        group_weights[[group_idx]],
        lambda,
        id_folds[group_ids]
      )
    }


    for (id_idx in 1:length(unique(ids))) {
      id_set <- ids == unique(ids)[id_idx]
      id_group  <- which(unique(groups) == unique(groups[id_set]))
      id_x[[id_idx]] <- x[id_set, ]
      id_params[[id_idx]] <- params[id_set, ]
      id_weights[[id_idx]] <- weights[id_set]
      id_scans <- droplevels(scans[id_set])
      if (length(unique(id_scans)) < k) {
        temp_k <- length(unique(id_scans))
      } else {
        temp_k <- k
      }
      scan_fold_vec <- rep(1:temp_k, ceiling(length(unique(id_scans)) / temp_k))
      scan_folds <- scan_fold_vec[
        sample(
          unique(id_scans),
          length(unique(id_scans)),
          replace = FALSE
        )
      ]
      scan_fold_vals <- scan_folds[id_scans]
      if (length(unique(id_scans)) == 1) {
        scan_fold_vec <- rep(1:k, ceiling(nrow(id_x[[id_idx]]) / k))
        scan_fold_vals <- scan_fold_vec[
          sample(
            1:nrow(id_x[[id_idx]]),
            nrow(id_x[[id_idx]]),
            replace = FALSE
          )
        ]
      }
      id_preds <- map(
        1:nrow(id_params[[id_idx]]),
        ~ (
          population_embedding$embedding_map(id_params[[id_idx]][.x, ]) +
          group_embeddings[[id_group]]$embedding_map(id_params[[id_idx]][.x, ]) +
          img_embeddings[[scans[id_set][.x]]]$embedding_map(id_params[[id_idx]][.x, ])
        )
      ) %>%
        reduce(rbind)
      id_embeddings[[id_idx]] <- fit_weighted_spline(
        id_x[[id_idx]] - id_preds,
        id_params[[id_idx]],
        id_weights[[id_idx]],
        lambda,
        scan_fold_vals
      )
    }

    for (img_idx in seq_along(unique(scans))) {
      img_set <- scans == unique(scans)[img_idx]
      img_group  <- which(unique(groups) == unique(groups[img_set]))
      img_id <- which(unique(ids) == unique(ids[img_set]))
      img_x[[img_idx]] <- x[img_set, ]
      img_params[[img_idx]] <- params[img_set, ]
      img_weights[[img_idx]] <- weights[img_set]
      obs_fold_vec <- rep(1:k, ceiling(sum(img_set) / k))
      obs_folds <- obs_fold_vec[
        sample(
          1:sum(img_set),
          sum(img_set),
          replace = FALSE
        )
      ]
      img_preds <- map(
        1:nrow(img_params[[img_idx]]),
        ~ (
          population_embedding$embedding_map(img_params[[img_idx]][.x, ]) +
          group_embeddings[[img_group]]$embedding_map(img_params[[img_idx]][.x, ]) +
          id_embeddings[[img_id]]$embedding_map(img_params[[img_idx]][.x, ])
        )
      ) %>%
        reduce(rbind)
      img_embeddings[[img_idx]] <- fit_weighted_spline(
        img_x[[img_idx]] - img_preds,
        img_params[[img_idx]],
        img_weights[[img_idx]],
        lambda,
        obs_folds
      )
      mean_x[[img_idx]] <- map(
        1:nrow(img_params[[img_idx]]),
        ~ (
          population_embedding$embedding_map(img_params[[img_idx]][.x, ]) +
          group_embeddings[[img_group]]$embedding_map(img_params[[img_idx]][.x, ]) +
          id_embeddings[[img_id]]$embedding_map(img_params[[img_idx]][.x, ]) +
          img_embeddings[[img_idx]]$embedding_map(img_params[[img_idx]][.x, ])
        )
      ) %>%
        reduce(rbind) %>%
        colMeans()
    }

    n <- n + 1
    # plot_id_pred(
    #   sample(unique(ids), 1),
    #   population_embedding,
    #   id_embeddings,
    #   params,
    #   ids,
    #   mean_x
    # )

    x_preds <- map(
      1:nrow(x),
      ~ (
          population_embedding$embedding_map(params[.x, ]) +
          group_embeddings[[groups[.x]]]$embedding_map(params[.x, ]) +
          id_embeddings[[ids[.x]]]$embedding_map(params[.x, ]) +
          img_embeddings[[scans[.x]]]$embedding_map(params[.x, ])
        )
    ) %>%
      reduce(rbind)

    mse <- map(1:nrow(x), ~ weights[.x] * dist_euclidean(x[.x, ], x_preds[.x, ])^2) %>%
      reduce(c) %>%
      mean()

    mse_ratio <- abs(mse - mse_old) / mse_old

    epsilon_hat <- mse_ratio
    print(
      paste0(
        "Maximum relative coefficient difference: ",
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
  for (img_idx in 1:length(unique(scans))) {
    scan_set <- scans == unique(scans)[img_idx]
    group_idx <- which(unique(groups) == unique(groups[scan_set]))
    id_idx <- which(unique(ids) == unique(ids[scan_set]))
    full_embeddings[[id_idx]] <- function(parameters) {
      population_embedding$embedding_map(parameters) +
        group_embeddings[[group_idx]]$embedding_map(parameters) +
        id_embeddings[[id_idx]]$embedding_map(parameters) +
        img_embeddings[[img_idx]]$embedding_map(parameters) -
        mean_x[[img_idx]]
    }
  }

  return(
    list(
      embeddings = full_embeddings,
      population_embedding = population_embedding,
      group_embeddings = group_embeddings,
      id_embeddings = id_embeddings,
      img_embeddings = img_embeddings
    )
  )
}
