library(doParallel)
library(foreach)
library(here)
library(pme)
library(idx2r)
library(progress)
library(Rfast)
library(tidyverse)

source(here("code/functions/fit_weighted_spline.R"))
source(here("code/functions/print_SSD.R"))
source(here("code/functions/fit_mnist_additive_model.R"))
source(here("code/functions/plot_id_pred.R"))

SSD_ratio_threshold <- 5
verbose <- TRUE


mnist_url <- "https://systemds.apache.org/assets/datasets/mnist/"
images_file <- "train-images-idx3-ubyte.gz"
labels_file <- "train-labels-idx1-ubyte.gz"
download.file(
  paste0(mnist_url, images_file),
  destfile = here(paste0("data/", images_file))
)
download.file(
  paste0(mnist_url, labels_file),
  destfile = here(paste0("data/", labels_file))
)

R.utils::gunzip(here(paste0("data/", images_file)))
R.utils::gunzip(here(paste0("data/", labels_file)))


mnist_train <- read_idx(here("data/train-images-idx3-ubyte"))
mnist_labels <- read_idx(here("data/train-labels-idx1-ubyte"))


mnist_train5 <- mnist_train[mnist_labels == 5, , ]
mnist_train5_ids <- seq_len(nrow(mnist_train5))


rotate_img <- function(x) {
  # rotate image matrix 90 degrees clockwise
  t(apply(x, 2, rev))
}


mnist_train5_proc <- list()
for (idx in mnist_train5_ids) {
  mnist_mat <- rotate_img(mnist_train5[idx, , ])

  digit_indices <- which(mnist_mat != 0, arr.ind = TRUE)

  digit_indices_range <- colMinsMaxs(digit_indices)
  digit_indices[, 1] <- (digit_indices[, 1] - digit_indices_range[1, 1]) /
    (digit_indices_range[2, 1] - digit_indices_range[1, 1])
  digit_indices[, 2] <- (digit_indices[, 2] - digit_indices_range[1, 2]) /
    (digit_indices_range[2, 2] - digit_indices_range[1, 2])
  digit_indices <- cbind(idx, digit_indices)
  mnist_train5_proc[[idx]] <- digit_indices
}

mnist_train5 <- reduce(mnist_train5_proc, rbind)

mnist_train3 <- mnist_train[mnist_labels == 3, , ]
mnist_train3_ids <- seq_len(nrow(mnist_train3))

mnist_train3_proc <- list()
for (idx in mnist_train3_ids) {
  mnist_mat <- rotate_img(mnist_train3[idx, , ])
  digit_indices <- which(mnist_mat != 0, arr.ind = TRUE)

  digit_indices_range <- colMinsMaxs(digit_indices)
  digit_indices[, 1] <- (digit_indices[, 1] - digit_indices_range[1, 1]) /
    (digit_indices_range[2, 1] - digit_indices_range[1, 1])
  digit_indices[, 2] <- (digit_indices[, 2] - digit_indices_range[1, 2]) /
    (digit_indices_range[2, 2] - digit_indices_range[1, 2])
  digit_indices <- cbind(idx, digit_indices)
  mnist_train3_proc[[idx]] <- digit_indices
}

mnist_train3 <- reduce(mnist_train3_proc, rbind)

mnist_train5 <- mnist_train5[mnist_train5[, 1] <= 200, ]

id_vals <- mnist_train5[, 1]
input <- mnist_train5[, -1]

init_pme <- pme(
  input[id_vals == 1, ],
  d = 1,
  lambda = exp(-5:10),
  print_plots = FALSE,
  verbose = TRUE
)

opt_run <- which.min(init_pme$MSD)

mnist_red <- list()
mnist_x <- list()
mnist_weights <- list()
mnist_id <- list()
for (id in unique(id_vals)) {
  print(id)
  mnist_red[[id]] <- hdmde(input[id_vals == id, ], 20, 0.05, 100)
  mnist_x[[id]] <- mnist_red[[id]]$mu
  mnist_weights[[id]] <- mnist_red[[id]]$theta_hat
  mnist_id[[id]] <- rep(id, nrow(mnist_x[[id]]))
}
mnist_x <- reduce(mnist_x, rbind)
mnist_weights <- reduce(mnist_weights, c)
mnist_id <- reduce(mnist_id, c)

nearest_cluster <- vector()
for (i in 1:nrow(mnist_x)) {
  distances <- map(
    1:nrow(init_pme$parameterization[[opt_run]]),
    ~ dist_euclidean(
      mnist_x[i, ],
      init_pme$embedding_map(init_pme$parameterization[[opt_run]][.x, ])
    )
  ) %>%
    reduce(c)
  nearest_cluster[i] <- which.min(distances)
}

init_params <- map(
  seq_len(nrow(mnist_x)),
  ~ projection_pme(
    mnist_x[.x, ],
    init_pme$embedding_map,
    init_pme$parameterization[[opt_run]][nearest_cluster[.x], ]
  )
) %>%
  reduce(c) %>%
  matrix(ncol = 1)

lambda <- c(0, exp(-15:5))
k <- 5
fold_vec <- rep(1:k, ceiling(length(unique(mnist_id)) / k))
id_folds <- fold_vec[sample(
  seq_along(unique(mnist_id)),
  length(unique(mnist_id)),
  replace = FALSE
)]
d <- 1

epsilon <- 0.05
max_iter <- 250


init_additive_mod <- fit_mnist_additive_model(
  x = mnist_x,
  params = init_params,
  weights = mnist_weights,
  lambda = lambda,
  k = 5,
  ids = mnist_id,
  epsilon = 0.01,
  max_iter = 250
)

params <- map(
  1:nrow(mnist_x),
  ~ projection_pme(
    mnist_x[.x, ],
    function(x) {
      init_additive_mod$population_embedding$embedding_map(x) +
        init_additive_mod$id_embeddings[[mnist_id[.x]]]$embedding_map(x)
    },
    init_params[.x, ]
  )
) %>%
  reduce(c) %>%
  matrix(nrow = nrow(mnist_x))

SSD <- map(
  1:nrow(mnist_x),
  ~ dist_euclidean(
    mnist_x[.x, ],
    init_additive_mod$population_embedding$embedding_map(params[.x, ]) +
      init_additive_mod$id_embeddings[[mnist_id[.x]]]$embedding_map(params[
        .x,
      ])
  )^2
) %>%
  reduce(c) %>%
  sum()

# embeddings <- init_additive_mod$embeddings

count <- 1
SSD_ratio <- 10 * epsilon

while (
  (SSD_ratio > epsilon) &
    (SSD_ratio <= SSD_ratio_threshold) &
    (count <= (max_iter - 1))
) {
  SSD_prev <- SSD
  additive_model_old <- init_additive_mod
  params_prev <- params

  additive_model <- fit_mnist_additive_model(
    mnist_x,
    params = params,
    weights = mnist_weights,
    lambda = lambda,
    k = 5,
    ids = mnist_id,
    epsilon = 0.01,
    max_iter = 250
  )

  params <- map(
    1:nrow(mnist_x),
    ~ projection_pme(
      mnist_x[.x, ],
      function(x) {
        additive_model$population_embedding$embedding_map(x) +
          additive_model$id_embeddings[[mnist_id[.x]]]$embedding_map(x)
      },
      params[.x, ]
    )
  ) %>%
    reduce(c) %>%
    matrix(nrow = nrow(mnist_x))

  SSD <- map(
    1:nrow(mnist_x),
    ~ dist_euclidean(
      mnist_x[.x, ],
      additive_model$population_embedding$embedding_map(params[.x, ]) +
        additive_model$id_embeddings[[mnist_id[.x]]]$embedding_map(params[.x, ])
    )^2
  ) %>%
    reduce(c) %>%
    sum()

  SSD_ratio <- abs(SSD - SSD_prev) / SSD_prev
  count <- count + 1

  if (verbose == TRUE) {
    print_SSD(SSD, SSD_ratio, count)
  }
}

nearest_cluster_full <- vector()
input_red <- unique(mnist_train5)
pb <- progress_bar$new(
  total = nrow(input_red),
  format = "[:bar] :percent eta: :eta",
  clear = FALSE
)
for (i in 1:nrow(input_red)) {
  pb$tick()
  distances <- map(
    1:nrow(params),
    ~ dist_euclidean(
      input_red[i, -1],
      additive_model$population_embedding$embedding_map(params[.x, ]) +
        additive_model$id_embeddings[[input_red[i, 1]]]$embedding_map(params[
          .x,
        ])
    )
  ) %>%
    reduce(c)
  nearest_cluster_full[i] <- which.min(distances)
}

input_df_names <- c("id", "x", "y")
mnist_train5_df <- data.frame(mnist_train5)
input_red_df <- data.frame(input_red)
names(mnist_train5_df) <- input_df_names
names(input_red_df) <- input_df_names
input_red_df$cluster <- nearest_cluster_full

mnist_train5_df <- full_join(
  mnist_train5_df,
  input_red_df,
  by = c("id", "x", "y")
)
mnist_train5_mat <- as.matrix(mnist_train5_df)

mnist_params <- map(
  1:nrow(mnist_train5_mat),
  ~ {
    print(.x)
    projection_pme(
      mnist_train5_mat[.x, 2:3],
      function(x) {
        additive_model$population_embedding$embedding_map(x) +
          additive_model$id_embeddings[[mnist_train5_mat[.x, 1]]]$embedding_map(
            x
          )
      },
      params[mnist_train5_mat[.x, 4], ]
    )
  }
) %>%
  reduce(rbind) %>%
  matrix(nrow = nrow(mnist_train5_mat))

hpme_msd <- map(
  1:nrow(mnist_train5_mat),
  ~ {
    print(.x)
    dist_euclidean(
      mnist_train5_mat[.x, 2:3],
      additive_model$population_embedding$embedding_map(mnist_params[.x, ]) +
        additive_model$id_embeddings[[mnist_train5_mat[
          .x,
          1
        ]]]$embedding_map(mnist_params[.x, ])
    )^2
  }
) %>%
  reduce(c) %>%
  mean()

mnist_train5_out <- list(
  data = mnist_train5,
  model = additive_model,
  params = params
)

saveRDS(mnist_train5_out, "mnist_train5_result.RDS")
