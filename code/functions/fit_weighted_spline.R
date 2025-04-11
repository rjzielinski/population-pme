fit_weighted_spline <- function(x, params, weights, smoothing_vals, folds, parallel = "sequential") {
  require(future, quietly = TRUE)
  require(pme, quietly = TRUE)
  require(purrr, quietly = TRUE)
  k <- length(unique(folds))
  coefs <- list()
  d <- ncol(params)
  smoothing_error <- vector()
  lambda <- smoothing_vals

  cores <- availableCores()

  for (smoothing_idx in seq_along(smoothing_vals)) {
    coefs[[smoothing_idx]] <- list()

    if (parallel == "sequential") {
      plan(sequential)
    } else if (parallel == "multisession") {
      plan(multisession, workers = min(cores, k))
    } else if (parallel == "multicore") {
      plan(multicore, workers = min(cores, k))
    }
    fold_output <- list()
    for (k_idx in seq_len(k)) {
      fold_output[[k_idx]] <- future({
        train_params <- params[!(folds == k_idx), ] %>%
          matrix(ncol = d)
        fold_params <- params[folds == k_idx, ] %>%
          matrix(ncol = d)
        train_param_mat <- cbind(1, train_params)
        train_x <- x[!(folds == k_idx), ]
        fold_x <- x[folds == k_idx, ]
        train_weights <- diag(weights[!(folds == k_idx)])
        train_E <- calcE(train_params, 4 - ncol(train_params))
        coefs[[smoothing_idx]][[k_idx]] <- solve_weighted_spline(
          train_E,
          train_weights,
          train_param_mat,
          train_x,
          lambda[smoothing_idx],
          ncol(train_params),
          ncol(train_x)
        )
        fold_embedding <- function(parameters) {
          as.vector(
            (t(coefs[[smoothing_idx]][[k_idx]][seq_len(nrow(train_x)), ]) %*% etaFunc(parameters, train_params, 4 - d)) +
              (t(coefs[[smoothing_idx]][[k_idx]][(nrow(train_x) + 1):(nrow(train_x) + d + 1), ]) %*% matrix(c(1, parameters), ncol = 1))
          )
        }
        error_val <- map(
          seq_len(nrow(fold_x)),
          ~ weights[folds == k_idx][.x] * dist_euclidean(
            fold_x[.x, ], fold_embedding(fold_params[.x, ])
          )^2
        ) %>%
          reduce(c) %>%
          mean()
        error_val
      }, seed = TRUE)
    }
    fold_errors <- map(fold_output, ~ value(.x)) %>%
      reduce(c)
    smoothing_error[smoothing_idx] <- mean(fold_errors)
    if (length(smoothing_error) >= 5) {
      if (!is.unsorted(smoothing_error[(length(smoothing_error) - 4):length(smoothing_error)])) {
        break
      }
    }
  }
  opt_run <- which.min(smoothing_error)

  E <- calcE(params, 4 - ncol(params))
  param_mat <- cbind(rep(1, nrow(params)), params)
  opt_coefs <- solve_weighted_spline(
    E,
    diag(weights),
    param_mat,
    x,
    smoothing_vals[opt_run],
    ncol(params),
    ncol(x)
  )
  embedding_map <- function(parameters) {
    as.vector(
      (t(opt_coefs[seq_len(nrow(x)), ]) %*% etaFunc(parameters, params, 4 - d)) +
        (t(opt_coefs[(nrow(x) + 1):(nrow(x) + d + 1), ]) %*% matrix(c(1, parameters), ncol = 1))
    )
  }
  return(
    list(
      embedding_map = embedding_map,
      coefs = opt_coefs,
      lambda = smoothing_vals[opt_run]
    )
  )
}
