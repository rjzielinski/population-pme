library(doParallel)
library(foreach)
# library(pme)
library(ggplot2)
library(ggimage)
library(imager)
library(progress)
library(Rfast)
library(tidyverse)

plot_img <- function(img_vec) {
  img_mat <- matrix(img_vec, nrow = 28)[, 28:1]
  image(img_mat)
}

mnist <- readRDS("~/Documents/brown/research/population-pme/data/mnist.RDS")
mnist_train5 <- mnist$train$images[mnist$train$labels == 5, ]
mnist_train5_ids <- 1:nrow(mnist_train5)

mnist_train <- mnist$train$images
mnist_train_ids <- 1:nrow(mnist_train)
mnist_train_labels <- mnist$train$labels

set.seed(58173)
mnist_train_idx <- sample(1:nrow(mnist_train), 2500)
mnist_train <- mnist_train[mnist_train_idx, ]
mnist_train_labels <- mnist_train_labels[mnist_train_idx]

mnist_pme <- pme(
  mnist_train,
  d = 3,
  initialization_type = "subsample",
  lambda = exp(-20:5),
  epsilon = 0.01,
  alpha = 0.01,
  max_iter = 1000,
  verbose = TRUE
)

params <- mnist_pme$parameterization[[which.min(mnist_pme$MSD)]]
est_centers <- purrr::map(
  1:nrow(params),
  ~ mnist_pme$embedding_map(params[.x, ])
) %>%
  reduce(rbind)

nearest_centers <- calc_nearest_x(mnist_train, est_centers)

proj_params <- purrr::map(
  1:nrow(mnist_train),
  ~ projection_pme(mnist_train[.x, ], f = mnist_pme$embedding_map, initial_guess = params[nearest_centers[.x] + 1, ])
) %>%
  reduce(rbind)


mnist_images <- list()
for (idx in 1:nrow(mnist_train)) {
  # mnist_mat <- matrix(mnist_train5[idx, ], nrow = 28)[, 28:1]
  # mnist_mat <- matrix(mnist_train5[idx, ], nrow = 28)[, 28:1]
  mnist_mat <- matrix(mnist_train[idx, ], nrow = 28)
  mnist_images[[idx]] <- as.cimg(mnist_mat)
}

for (img_val in 1:length(mnist_images)) {
  image_path <- paste0("~/Documents/brown/research/population-pme/data/mnist_images/img_", str_pad(img_val, 4, side = "left", pad = "0"), ".png")
  imager::save.image(mnist_images[[img_val]], image_path)
}

img <- list.files(
  path="~/Documents/brown/research/population-pme/data/mnist_images",
  pattern="png",
  full.names=TRUE
)

img_df <- data.frame(
  x = proj_params[, 1],
  y = proj_params[, 2],
  z = proj_params[, 3],
  image = img
)

img_df_test <- img_df[1:100, ]

ggplot(img_df_test, aes(x = x, y = y, color = z)) + geom_image(aes(image=image), size=.05)


# mnist_train5_mat <- mnist_train5[1:100, ]

mnist_train5_pme <- pme(mnist_train5, d = 2, verbose = TRUE)
params <- mnist_train5_pme$parameterization[[which.min(mnist_train5_pme$MSD)]]
est_centers <- purrr::map(
  1:nrow(params),
  ~ mnist_train5_pme$embedding_map(params[.x, ])
) %>%
  reduce(rbind)

nearest_centers <- calc_nearest_x(mnist_train5_mat, est_centers)

proj_params <- purrr::map(
  1:nrow(mnist_train5_mat),
  ~ projection_pme(mnist_train5_mat[.x, ], f = mnist_train5_pme$embedding_map, initial_guess = params[nearest_centers[.x] + 1, ])
) %>%
  reduce(rbind)

# mnist_train_images <- mnist_images[1:100]

for (img_val in 1:length(mnist_train_images)) {
  image_path <- paste0("~/Documents/brown/research/population-pme/data/mnist_images/img_", str_pad(img_val, 4, side = "left", pad = "0"), ".png")
  imager::save.image(mnist_train_images[[img_val]], image_path)
}

img <- list.files(
  path="~/Documents/brown/research/population-pme/data/mnist_images",
  pattern="png",
  full.names=TRUE
)

img_df <- data.frame(
  x = proj_params[, 1],
  y = proj_params[, 2],
  image = img
)

img_df_test <- img_df[1:20, ]

ggplot(img_df_test, aes(x, y)) + geom_image(aes(image=image), size=.05)
