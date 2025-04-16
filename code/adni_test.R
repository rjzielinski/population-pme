library(doParallel)
library(doFuture)
library(foreach)
library(listenv)
library(pme)
library(progress)
library(Rfast)
library(tidyverse)

options(future.globals.maxSize = 32 * 1024^3)
cores <- detectCores() - 1

source("code/functions/fit_weighted_spline.R")
source("code/functions/print_SSD.R")
source("code/functions/fit_adni_additive_model.R")
source("code/functions/plot_id_pred.R")
source("code/functions/calculate_pme_reconstructions.R")

SSD_ratio_threshold <- 5
verbose <- TRUE

lhipp_surface <- read_csv("data_old/adni_fsl_lhipp_surface.csv")
adni_info <- read_csv("data_old/adni_fsl_info.csv")

lhipp_surface <- lhipp_surface |>
  mutate(
    date = decimal_date(date),
    Group = gsub("EMCI", "MCI", Group),
    Group = gsub("LMCI", "MCI", Group),
  )

lhipp_surface_centers <- lhipp_surface |>
  group_by(subid, date, scan_id) |>
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x - mean_x)),
    max_y = max(abs(y - mean_y)),
    max_z = max(abs(z - mean_z))
  )

lhipp_bl <- lhipp_surface |>
  group_by(subid) |>
  arrange(date) |>
  summarize(
    date_bl = first(date),
    max_date = max(date)
  ) |>
  mutate(duration = max_date - date_bl)

lhipp_surface <- lhipp_surface |>
  full_join(lhipp_surface_centers, by = c("subid", "date", "scan_id")) |>
  full_join(lhipp_bl, by = "subid") |>
  mutate(
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    time_from_bl = date - date_bl,
    partition = ifelse(x > 0, 1, 2)
  )

lhipp_surface_inputs <- lhipp_surface |>
  select(date, x, y, z) |>
  as.matrix()

lhipp_groups <- as.factor(lhipp_surface$Group)
lhipp_ids <- as.factor(lhipp_surface$subid)
lhipp_scans <- as.factor(lhipp_surface$image_id)
lhipp_partition <- lhipp_surface$partition

lhipp_init_data <- list()
lhipp_init_data[[1]] <- lhipp_surface_inputs[
  (lhipp_groups == lhipp_groups[1]) &
    (lhipp_ids == lhipp_ids[1]) &
    (lhipp_scans == lhipp_scans[1]) &
    (lhipp_partition == 1),
]
lhipp_init_data[[2]] <- lhipp_surface_inputs[
  (lhipp_groups == lhipp_groups[1]) &
    (lhipp_ids == lhipp_ids[1]) &
    (lhipp_scans == lhipp_scans[1]) &
    (lhipp_partition == 2),
]

plan(multisession, workers = 2)
init_pme_wrapper <- list()
for (i in seq_len(2)) {
  init_pme_wrapper[[i]] <- future(
    {
      pme(
        lhipp_init_data[[i]][, -1],
        d = 2,
        print_plots = FALSE,
        verbose = FALSE
      )
    },
    seed = TRUE
  )
}

lhipp_init_pme <- map(
  init_pme_wrapper,
  ~ value(.x)
)

plan(multisession, workers = cores - 1)
lhipp_reduced <- list()
scan_list <- unique(lhipp_scans)
for (scan_idx in seq_along(scan_list)) {
  lhipp_reduced[[scan_idx]] <- future(
    {
      scan_indices <- lhipp_scans == scan_list[scan_idx]
      group_val <- unique(lhipp_groups[scan_indices])
      id_val <- unique(lhipp_ids[scan_indices])

      scan_data_pt1 <- lhipp_surface_inputs[
        (lhipp_groups == group_val) &
          (lhipp_ids == id_val) &
          (lhipp_scans == scan_list[scan_idx]) &
          (lhipp_partition == 1),
      ]
      scan_time <- scan_data_pt1[1, 1]
      scan_data_pt1 <- scan_data_pt1[, -1]

      scan_data_pt2 <- lhipp_surface_inputs[
        (lhipp_groups == group_val) &
          (lhipp_ids == id_val) &
          (lhipp_scans == scan_list[scan_idx]) &
          (lhipp_partition == 2), -1
      ]

      scan_red_pt1 <- hdmde(scan_data_pt1, 10, 0.05, 100)
      scan_red_pt2 <- hdmde(scan_data_pt2, 10, 0.05, 100)

      scan_x_pt1 <- scan_red_pt1$mu
      scan_x_pt2 <- scan_red_pt2$mu
      scan_weights_pt1 <- scan_red_pt1$theta_hat
      scan_weights_pt2 <- scan_red_pt2$theta_hat
      n_pt1 <- nrow(scan_x_pt1)
      n_pt2 <- nrow(scan_x_pt2)

      list(
        lhipp_x_pt1 = scan_x_pt1,
        lhipp_x_pt2 = scan_x_pt2,
        lhipp_weights_pt1 = scan_weights_pt1,
        lhipp_weights_pt2 = scan_weights_pt2,
        group_vec_pt1 = rep(group_val, n_pt1),
        group_vec_pt2 = rep(group_val, n_pt2),
        id_vec_pt1 = rep(id_val, n_pt1),
        id_vec_pt2 = rep(id_val, n_pt2),
        scan_vec_pt1 = rep(scan_list[scan_idx], n_pt1),
        scan_vec_pt2 = rep(scan_list[scan_idx], n_pt2),
        lhipp_time_pt1 = rep(scan_time, n_pt1),
        lhipp_time_pt2 = rep(scan_time, n_pt2)
      )
    },
    seed = TRUE
  )
}

group_vecs <- list(list(), list())
id_vecs <- list(list(), list())
img_vecs <- list(list(), list())
lhipp_times <- list(list(), list())

lhipp_x <- list(list(), list())
lhipp_weights <- list(list(), list())



for (scan_idx in seq_along(lhipp_reduced)) {
  lhipp_x[[1]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$lhipp_x_pt1
  lhipp_x[[2]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$lhipp_x_pt2
  lhipp_weights[[1]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$lhipp_weights_pt1
  lhipp_weights[[2]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$lhipp_weights_pt2
  group_vecs[[1]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$group_vec_pt1
  group_vecs[[2]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$group_vec_pt2
  id_vecs[[1]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$id_vec_pt1
  id_vecs[[2]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$id_vec_pt2
  img_vecs[[1]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$scan_vec_pt1
  img_vecs[[2]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$scan_vec_pt2
  lhipp_times[[1]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$lhipp_time_pt1
  lhipp_times[[2]][[scan_idx]] <- value(lhipp_reduced[[scan_idx]])$lhipp_time_pt2
}

lhipp_x[[1]] <- reduce(lhipp_x[[1]], rbind)
lhipp_x[[2]] <- reduce(lhipp_x[[2]], rbind)
lhipp_weights[[1]] <- reduce(lhipp_weights[[1]], c)
lhipp_weights[[2]] <- reduce(lhipp_weights[[2]], c)
group_vecs[[1]] <- reduce(group_vecs[[1]], c)
group_vecs[[2]] <- reduce(group_vecs[[2]], c)
id_vecs[[1]] <- reduce(id_vecs[[1]], c)
id_vecs[[2]] <- reduce(id_vecs[[2]], c)
img_vecs[[1]] <- reduce(img_vecs[[1]], c)
img_vecs[[2]] <- reduce(img_vecs[[2]], c)
lhipp_times[[1]] <- reduce(lhipp_times[[1]], c)
lhipp_times[[2]] <- reduce(lhipp_times[[2]], c)

init_params <- list()
init_params[[1]] <- calculate_pme_reconstructions(
  lhipp_init_pme[[1]],
  # lhipp_surface_inputs[lhipp_partition == 1, -1]
  lhipp_x[[1]]
)$projections
init_params[[2]] <- calculate_pme_reconstructions(
  lhipp_init_pme[[2]],
  lhipp_x[[2]]
)$projections

weights <- list()
weights[[1]] <- diag(lhipp_weights[[1]])
weights[[2]] <- diag(lhipp_weights[[2]])

init_additive_mod_pt1 <- fit_adni_additive_model(
  x = lhipp_x[[1]],
  params = init_params[[1]],
  weights = lhipp_weights[[1]],
  lambda = exp(-20:5),
  k = 5,
  groups = group_vecs[[1]],
  ids = id_vecs[[1]],
  scans = img_vecs[[1]],
  times = lhipp_times[[1]],
  epsilon = 0.01,
  max_iter = 250,
  cores = cores
)
init_additive_mod_pt2 <- fit_adni_additive_model(
  x = lhipp_x[[2]],
  params = init_params[[2]],
  weights = lhipp_weights[[2]],
  lambda = exp(-20:5),
  k = 5,
  groups = group_vecs[[2]],
  ids = id_vecs[[2]],
  scans = img_vecs[[2]],
  times = lhipp_times[[2]],
  epsilon = 0.01,
  max_iter = 250,
  cores = cores - 1
)

lambda <- exp(-20:5)
k <- 5
id_fold_vec <- rep(1:k, ceiling(length(unique(lhipp_ids)) / k))
id_folds <- id_fold_vec[sample(
  seq_along(unique(lhipp_ids)),
  length(unique(lhipp_ids)),
  replace = FALSE
)]

img_fold_vec <- rep(1:k, ceiling(length(unique(lhipp_scans)) / k))
img_folds <- img_fold_vec[sample(
  seq_along(unique(lhipp_scans)),
  length(unique(lhipp_scans)),
  replace = FALSE
)]

epsilon <- 0.05
max_iter <- 250

params <- list()
params[[1]] <- map(
  seq_len(nrow(lhipp_x[[1]])),
  ~ projection_pme(
    lhipp_x[[1]][.x, ],
    function(x) {
      init_additive_mod_pt1$population_embedding$embedding_map(x) +
        init_additive_mod_pt1$group_embeddings[[group_vecs[[1]][.x]]]$embedding_map(x) +
        init_additive_mod_pt1$id_embeddings[[id_vecs[[1]][.x]]]$embedding_map(x) +
        init_additive_mod_pt1$img_embeddings[[img_vecs[[1]][.x]]]$embedding_map(x)
    },
    init_params[[1]][.x, ]
  )
) %>%
  reduce(rbind)
params[[2]] <- map(
  seq_len(nrow(lhipp_x[[2]])),
  ~ projection_pme(
    lhipp_x[[2]][.x, ],
    function(x) {
      init_additive_mod_pt2$population_embedding$embedding_map(x) +
        init_additive_mod_pt2$group_embeddings[[group_vecs[[2]][.x]]]$embedding_map(x) +
        init_additive_mod_pt2$id_embeddings[[id_vecs[[2]][.x]]]$embedding_map(x) +
        init_additive_mod_pt2$img_embeddings[[img_vecs[[2]][.x]]]$embedding_map(x)
    },
    init_params[[2]][.x, ]
  )
) %>%
  reduce(rbind)

SSD_pt1 <- map(
  seq_len(nrow(lhipp_x[[1]])),
  ~ dist_euclidean(
    lhipp_x[[1]][.x, ],
    init_additive_mod_pt1$population_embedding$embedding_map(params[[1]][.x, ]) +
      init_additive_mod_pt1$group_embeddings[[group_vecs[[1]][.x]]]$embedding_map(params[[1]][.x, ]) +
      init_additive_mod_pt1$id_embeddings[[id_vecs[[1]][.x]]]$embedding_map(params[[1]][.x, ]) +
      init_additive_mod_pt1$img_embeddings[[img_vecs[[1]][.x]]]$embedding_map(params[[1]][.x, ])
  )^2
) %>%
  reduce(c) %>%
  sum()
SSD_pt2 <- map(
  seq_len(nrow(lhipp_x[[2]])),
  ~ dist_euclidean(
    lhipp_x[[2]][.x, ],
    init_additive_mod_pt2$population_embedding$embedding_map(params[[2]][.x, ]) +
      init_additive_mod_pt2$group_embeddings[[group_vecs[[2]][.x]]]$embedding_map(params[[2]][.x, ]) +
      init_additive_mod_pt2$id_embeddings[[id_vecs[[2]][.x]]]$embedding_map(params[[2]][.x, ]) +
      init_additive_mod_pt2$img_embeddings[[img_vecs[[2]][.x]]]$embedding_map(params[[2]][.x, ])
  )^2
) %>%
  reduce(c) %>%
  sum()

SSD <- SSD_pt1 + SSD_pt2

# embeddings <- init_additive_mod$embeddings

count <- 1
SSD_ratio <- 10 * epsilon

additive_model <- list()
additive_model$partition1 <- init_additive_mod_pt1
additive_model$partition2 <- init_additive_mod_pt2
additive_model$SSD <- c(SSD)

while ((SSD_ratio > epsilon) & (SSD_ratio <= SSD_ratio_threshold) & (count <= (max_iter - 1))) {
  SSD_prev <- SSD
  additive_model_old <- additive_model
  params_prev <- params

  additive_model$partition1 <- fit_adni_additive_model(
    x = lhipp_x[[1]],
    params = params[[1]],
    weights = lhipp_weights[[1]],
    lambda = exp(-20:5),
    k = 5,
    groups = group_vecs[[1]],
    ids = id_vecs[[1]],
    scans = img_vecs[[1]],
    times = lhipp_times[[1]],
    epsilon = 0.01,
    max_iter = 250,
    cores = cores - 1
  )
  additive_model$partition2 <- fit_adni_additive_model(
    x = lhipp_x[[2]],
    params = params[[2]],
    weights = lhipp_weights[[2]],
    lambda = exp(-20:5),
    k = 5,
    groups = group_vecs[[2]],
    ids = id_vecs[[2]],
    scans = img_vecs[[2]],
    times = lhipp_times[[2]],
    epsilon = 0.01,
    max_iter = 250,
    cores = cores - 1
  )

  params[[1]] <- map(
    seq_len(nrow(lhipp_x[[1]])),
    ~ projection_pme(
      lhipp_x[[1]][.x, ],
      function(x) {
        additive_model$partition1$population_embedding$embedding_map(x) +
          additive_model$partition1$group_embeddings[[group_vecs[[1]][.x]]]$embedding_map(x) +
          additive_model$partition1$id_embeddings[[id_vecs[[1]][.x]]]$embedding_map(x) +
          additive_model$partition1$img_embeddings[[img_vecs[[1]][.x]]]$embedding_map(x)
      },
      params_prev[[1]][.x, ]
    )
  ) %>%
    reduce(rbind)
  params[[2]] <- map(
    seq_len(nrow(lhipp_x[[2]])),
    ~ projection_pme(
      lhipp_x[[2]][.x, ],
      function(x) {
        additive_model$partition2$population_embedding$embedding_map(x) +
          additive_model$partition2$group_embeddings[[group_vecs[[2]][.x]]]$embedding_map(x) +
          additive_model$partition2$id_embeddings[[id_vecs[[2]][.x]]]$embedding_map(x) +
          additive_model$partition2$img_embeddings[[img_vecs[[2]][.x]]]$embedding_map(x)
      },
      params_prev[[2]][.x, ]
    )
  ) %>%
    reduce(rbind)



  SSD_pt1 <- map(
    seq_len(nrow(lhipp_x[[1]])),
    ~ dist_euclidean(
      lhipp_x[[1]][.x, ],
      additive_model$partition1$population_embedding$embedding_map(params[[1]][.x, ]) +
        additive_model$partition1$group_embeddings[[group_vecs[[1]][.x]]]$embedding_map(params[[1]][.x, ]) +
        additive_model$partition1$id_embeddings[[id_vecs[[1]][.x]]]$embedding_map(params[[1]][.x, ]) +
        additive_model$partition1$img_embeddings[[img_vecs[[1]][.x]]]$embedding_map(params[[1]][.x, ])
    )^2
  ) %>%
    reduce(c) %>%
    sum()
  SSD_pt2 <- map(
    seq_len(nrow(lhipp_x[[2]])),
    ~ dist_euclidean(
      lhipp_x[[2]][.x, ],
      additive_model$partition2$population_embedding$embedding_map(params[[2]][.x, ]) +
        additive_model$partition2$group_embeddings[[group_vecs[[2]][.x]]]$embedding_map(params[[2]][.x, ]) +
        additive_model$partition2$id_embeddings[[id_vecs[[2]][.x]]]$embedding_map(params[[2]][.x, ]) +
        additive_model$partition2$img_embeddings[[img_vecs[[2]][.x]]]$embedding_map(params[[2]][.x, ])
    )^2
  ) %>%
    reduce(c) %>%
    sum()

  SSD <- SSD_pt1 + SSD_pt2

  SSD_ratio <- abs(SSD - SSD_prev) / SSD_prev
  count <- count + 1
  additive_model$SSD[count] <- SSD

  if (verbose == TRUE) {
    print_SSD(SSD, SSD_ratio, count)
  }
}

nearest_cluster_full <- vector()
nearest_clusters <- matrix(
  nrow = nrow(lhipp_surface_inputs),
  ncol = ncol(lhipp_x[[1]])
)
nearest_params <- matrix(
  nrow = nrow(lhipp_surface_inputs),
  ncol = ncol(params[[1]])
)
pb <- progress_bar$new(
  total = nrow(lhipp_surface_inputs),
  format = "[:bar] :percent eta: :eta",
  clear = FALSE
)

for (i in seq_len(nrow(lhipp_surface_inputs))) {
  pb$tick()
  temp_clusters <- lhipp_x[[lhipp_partition[i]]][
    (group_vecs[[lhipp_partition[i]]] == as.numeric(lhipp_groups[i])) &
      (id_vecs[[lhipp_partition[i]]] == as.numeric(lhipp_ids[i])) &
      (img_vecs[[lhipp_partition[i]]] == as.numeric(lhipp_scans[i])),
  ]
  temp_params <- params[[lhipp_partition[i]]][
    (group_vecs[[lhipp_partition[i]]] == as.numeric(lhipp_groups[i])) &
      (id_vecs[[lhipp_partition[i]]] == as.numeric(lhipp_ids[i])) &
      (img_vecs[[lhipp_partition[i]]] == as.numeric(lhipp_scans[i])),
  ]
  distances <- map(
    seq_len(nrow(temp_clusters)),
    ~ dist_euclidean(
      lhipp_surface_inputs[i, -1],
      temp_clusters[.x, ]
    )
  ) |>
    reduce(c)
  nearest_cluster_full[i] <- which.min(distances)
  nearest_clusters[i, ] <- unlist(temp_clusters[which.min(distances), ])
  nearest_params[i, ] <- unlist(temp_params[which.min(distances), ])
}

lhipp_params <- map(
  seq_len(nrow(lhipp_surface_inputs)),
  ~ projection_pme(
    unlist(lhipp_surface_inputs[.x, -1]),
    function(x) {
      additive_model[[lhipp_partition[.x]]]$population_embedding$embedding_map(x) +
        additive_model[[lhipp_partition[.x]]]$group_embeddings[[as.numeric(lhipp_groups[.x])]]$embedding_map(x) +
        additive_model[[lhipp_partition[.x]]]$id_embeddings[[as.numeric(lhipp_ids[.x])]]$embedding_map(x) +
        additive_model[[lhipp_partition[.x]]]$img_embeddings[[as.numeric(lhipp_scans[.x])]]$embedding_map(x)
    },
    nearest_params[.x, ]
  ),
  .progress = TRUE
) |>
  reduce(rbind)

lhipp_projections <- map(
  seq_len(nrow(lhipp_surface_inputs)),
  ~ additive_model[[lhipp_partition[.x]]]$population_embedding$embedding_map(lhipp_params[.x, ]) +
    additive_model[[lhipp_partition[.x]]]$group_embeddings[[as.numeric(lhipp_groups[.x])]]$embedding_map(lhipp_params[.x, ]) +
    additive_model[[lhipp_partition[.x]]]$id_embeddings[[as.numeric(lhipp_ids[.x])]]$embedding_map(lhipp_params[.x, ]) +
    additive_model[[lhipp_partition[.x]]]$img_embeddings[[as.numeric(lhipp_scans[.x])]]$embedding_map(lhipp_params[.x, ]),
  .progress = TRUE
) |>
  reduce(rbind)

lhipp_projections_pop <- map(
  seq_len(nrow(lhipp_surface_inputs)),
  ~ additive_model[[lhipp_partition[.x]]]$population_embedding$embedding_map(lhipp_params[.x, ]),
  .progress = TRUE
) |>
  reduce(rbind)

hpme_msd <- map(
  seq_len(nrow(lhipp_surface_inputs)),
  ~ {
    dist_euclidean(
      lhipp_surface_inputs[.x, -1],
      lhipp_projections[.x, ]
    )^2
  }
) %>%
  reduce(c) %>%
  mean()

lhipp_test_out <- list(
  data = lhipp_surface,
  data_red = list(
    x = lhipp_x,
    weights = lhipp_weights,
    groups = group_vecs,
    ids = id_vecs,
    imgs = img_vecs
  ),
  model = additive_model,
  params = lhipp_params,
  msd = hpme_msd
)
saveRDS(lhipp_test_out, "output/test_lhipp_additive_model.RDS")
