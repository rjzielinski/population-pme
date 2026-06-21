#' Simulate Data
#'
#' Generates simulated data from given underlying manifold.
#'
#' @param duration The simulated study's follow-up duration
#' @param interval The scheduled interval of study visits
#' @param case Specifies the manifold used to generate the data
#' @param obs_noise The standard deviation of the independent Gaussian observation-level noise
#' @param amplitude_noise The standard deviation of Gaussian noise applied to the amplitude of the manifold at each time point
#' @param period_noise The standard deviation of Gaussian noise applied to the period of the manifold at each time point
#' @param time_trend Specifies systematic longitudinal changes in the manifold amplitude
#' @param time_change Scales the longitudinal amplitude changes
#' @param N The number of observations generated at each time point
#'
simulate_data <- function(
  duration,
  interval,
  case,
  obs_noise,
  amplitude_noise,
  period_noise,
  time_trend,
  time_change,
  visit_noise = 0.1,
  N = 1000
) {
  # instead of using set interval between visits consider using a parameter
  # for the expected number of follow up visits within the study
  # then draw the observed number of visits from a poisson distribution
  # visit times will then be uniformly distributed between 0 and the duration

  # simulate_data <- function(duration, exp_follow_up, case, obs_noise, amplitude_noise, period_noise, time_trend, time_change, N = 1000) {

  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  require(tibble, quietly = TRUE, warn.conflicts = FALSE)

  # create list of functions used to define manifolds of interest
  manifolds <- list(
    # case 1: sine curve with limited domain
    function(tau, amplitude, period) {
      c(tau, amplitude[1] * sin(period[1] * tau + pi / 2))
    },
    # case 2: sine curve with extended domain
    function(tau, amplitude, period) {
      c(tau, amplitude[1] * sin(period[1] * tau))
    },
    # case 3: partial circle
    function(tau, amplitude, period) {
      c(
        amplitude[1] * cos(period[1] * tau),
        amplitude[2] * sin(period[2] * tau)
      )
    },
    # case 4: polynomial curve
    function(tau, amplitude, period) {
      # note, the "period" argument name is maintained for consistency
      c(
        tau,
        (tau * amplitude[1] + period[1])^2,
        (tau * amplitude[2] + period[2])^3
      )
    },
    # case 5: corkscrew-like shape
    function(tau, amplitude, period) {
      c(
        tau,
        amplitude[1] * cos(period[1] * tau),
        amplitude[2] * sin(period[2] * tau)
      )
    },
    # case 6: bowl-shaped structure
    function(tau, amplitude, period) {
      # note that tau is 2-dimensional in this and following cases
      c(
        (period[1] / 2) * tau[1],
        (period[2] / 2) * tau[2],
        amplitude[1] *
          amplitude[2] *
          norm(
            matrix(period[1:d] * tau),
            type = "F"
          )
      )
    },
    # case 7: swiss roll
    function(tau, amplitude, period) {
      c(
        (period[1] * amplitude[1] * tau[1]) * cos(amplitude[1] * tau[1]),
        (period[2] * amplitude[2] * tau[1]) * sin(amplitude[2] * tau[1]),
        tau[2]
      )
    },
    # case 8: sphere
    function(tau, amplitude, period) {
      c(
        amplitude[1] * sin(period[1] * tau[1]) * cos(period[2] * tau[2]),
        amplitude[1] * sin(period[1] * tau[1]) * sin(period[2] * tau[2]),
        amplitude[1] * cos(period[1] * tau[1])
      )
    },
    # case 9: sphere
    function(tau, amplitude, period) {
      c(
        amplitude[1] * sin(period[1] * tau[1]) * cos(period[2] * tau[2]),
        amplitude[1] * sin(period[1] * tau[1]) * sin(period[2] * tau[2]),
        amplitude[1] * cos(period[1] * tau[1])
      )
    }
  )

  # specify dimensionality for all manifold cases
  D <- if (case == 1) {
    2
  } else if (case == 2) {
    2
  } else if (case == 3) {
    2
  } else if (case == 4) {
    3
  } else if (case == 5) {
    3
  } else if (case == 6) {
    3
  } else if (case == 7) {
    3
  } else if (case == 8) {
    3
  } else if (case == 9) {
    3
  }

  d <- if (case == 1) {
    1
  } else if (case == 2) {
    1
  } else if (case == 3) {
    1
  } else if (case == 4) {
    1
  } else if (case == 5) {
    1
  } else if (case == 6) {
    2
  } else if (case == 7) {
    2
  } else if (case == 8) {
    2
  } else if (case == 9) {
    2
  }

  # this assumes that study visits are evenly spaced.
  # what if this is not the case?
  # noise could be added to the time from baseline
  time_points <- seq(from = 0, to = duration, by = interval)

  # TODO: either incorporate or remove before publication
  #
  # case if using exp_follow_up
  # n_follow_up <- rpois(1, exp_follow_up)
  # time_points <- c(0, runif(n_follow_up, min = 0, max = duration))
  time_noise <- rnorm(n = length(time_points) - 1, mean = 0, sd = visit_noise)
  time_noise <- c(0, time_noise)
  # ensure
  time_points <- map(
    seq_along(time_points),
    ~ max(time_points[.x] + time_noise[.x], rnorm(1, 1e-2, 1e-3))
  ) |>
    reduce(c) |>
    sort()
  time_points[1] <- 0

  # specify systematic changes in amplitude over time
  time_adjustments <- if (time_trend == "constant") {
    0 * time_points
  } else if (time_trend == "linear") {
    time_points
  } else if (time_trend == "quadratic") {
    time_points^2
  }

  # scale time adjustments to be proportion of 1
  if (max(time_adjustments) == 0) {
    time_adjustments_scaled <- time_adjustments * time_change * 0
  } else {
    time_adjustments_scaled <- (time_adjustments / max(time_adjustments)) *
      time_change
  }
  # subtract from initial amplitude to replicate effects of atrophy
  amplitude_values <- replicate(2, 1 - time_adjustments_scaled)
  # for now assume that period does not change between visits
  period_values <- replicate(2, rep(1, length(time_points)))

  observed_amplitude_values <- list()
  observed_period_values <- list()

  sim_matrix <- matrix(NA, nrow = 1, ncol = 1 + d + (2 * D))
  manifold <- manifolds[[case]]

  for (time_idx in seq_along(time_points)) {
    # matrix of low-dimensional parameters
    r <- matrix(NA, nrow = N, ncol = d)

    # set domains for each of the manifolds
    if (case == 1) {
      r[, 1] <- rnorm(N, mean = 0, sd = 1)
    } else if (case == 2) {
      r[, 1] <- runif(N, min = -3 * pi, max = 3 * pi)
    } else if (case == 3) {
      r[, 1] <- runif(N, min = -0.9 * pi, max = 0 * pi)
    } else if (case == 4) {
      r[, 1] <- runif(N, min = -1, max = 1)
    } else if (case == 5) {
      r[, 1] <- runif(N, min = 0, max = 3 * pi)
    } else if (case == 6) {
      r[, 1] <- runif(N, min = -1, max = 1)
      r[, 2] <- runif(N, min = -1, max = 1)
    } else if (case == 7) {
      r[, 1] <- runif(N, min = (0.2 * pi), max = (2.9 * pi))
      r[, 2] <- runif(N, min = -1, max = 1)
    } else if (case == 8) {
      r[, 1] <- runif(N, min = 0, max = pi)
      r[, 2] <- runif(N, min = 0, max = 2 * pi)
    } else if (case == 9) {
      r[, 1] <- runif(N, min = 0, max = pi)
      r[, 2] <- runif(N, min = 0, max = 2 * pi)
    }

    # image-specific amplitude and period noise indicates measurement error
    # and variability in image processing/segmentation
    amplitude_noise_vals <- rnorm(d, mean = 0, sd = amplitude_noise)
    period_noise_vals <- rnorm(d, mean = 0, sd = period_noise)

    observed_amplitude_values[[time_idx]] <- map(
      seq_along(amplitude_values[time_idx, ]),
      ~ max(
        amplitude_values[time_idx, ][.x] +
          (amplitude_noise_vals * (1 - time_adjustments_scaled[time_idx])),
        rnorm(1, 1e-2, 1e-3)
      )
    ) |>
      reduce(c)

    observed_period_values[[time_idx]] <- map(
      seq_along(period_values[time_idx, ]),
      ~ max(
        period_values[time_idx, ][.x] +
          period_noise_vals,
        rnorm(1, 1e-2, 1e-3)
      )
    ) |>
      reduce(c)

    # for now assume that voxel-level errors are spatially-independent
    # to incorporate spatial dependencies, allow sd to vary by
    # parameter values
    # obs_noise_vals <- rnorm(N * D, mean = 0, sd = obs_noise) |>
    #   matrix(nrow = N, ncol = D)

    # test assumption that the level of voxel-level measurement error varies
    # spatially by simulating GP to determine variance of random noise
    gp_sigma <- calculate_sigma(r, d, obs_noise / 4, bandwidth = 1)
    obs_noise_adjustment <- MASS::mvrnorm(
      n = 1,
      mu = rep(0, N),
      Sigma = gp_sigma
    )
    obs_noise_adjusted <- (obs_noise) * (1 + obs_noise_adjustment)
    obs_noise_adjusted <- rep(obs_noise_adjusted, each = D)
    obs_noise_vals <- rnorm(N * D, mean = 0, sd = obs_noise_adjusted) |>
      matrix(nrow = N, ncol = D, byrow = TRUE)

    if (case == 9) {
      obs_noise_vals <- abs(obs_noise_vals)
    }

    # use true underlying amplitude and period values to generate
    # denoised observations
    X_true <- map(
      seq_len(nrow(r)),
      ~ manifold(
        r[.x, ],
        amplitude_values[time_idx, ],
        period_values[time_idx, ]
      )
    ) |>
      unlist() |>
      matrix(ncol = D, byrow = TRUE)

    # incorporate noise for manifold amplitude and period, and
    # observation-level noise for observations used to fit manifold learning
    X <- map(
      seq_len(nrow(r)),
      ~ manifold(
        r[.x, ],
        observed_amplitude_values[[time_idx]],
        observed_period_values[[time_idx]]
      )
    ) |>
      unlist() |>
      matrix(ncol = D, byrow = TRUE)

    X <- X + obs_noise_vals

    time_mat <- matrix(time_points[time_idx], nrow = nrow(r), ncol = 1)

    sim_data_mat <- cbind(time_mat, r, X, X_true)
    sim_matrix <- rbind(sim_matrix, sim_data_mat)
  }

  # remove initial row of NA values
  sim_matrix <- sim_matrix[-1, ]
  sim_df <- as_tibble(sim_matrix, .name_repair = "minimal")
  names(sim_df) <- c(
    "time",
    paste0("r", 1:d),
    paste0("X", 1:D),
    paste0("X_true", 1:D)
  )

  simulation_out <- list(
    df = sim_df,
    amplitude_values = observed_amplitude_values,
    period_values = observed_period_values
  )

  return(simulation_out)
}

calculate_sigma <- function(r, d, sd, bandwidth) {
  # creates a covariance matrix using the squared exponential kernel

  if (d == 1) {
    sigma_mat <- matrix(NA, nrow = nrow(r), ncol = nrow(r))
    for (i in seq_len(nrow(r))) {
      for (j in seq_len(nrow(r))) {
        dist <- r[i, 1] - r[j, 1]
        sigma_mat[i, j] <- sd^2 * exp(-(dist^2) / (2 * bandwidth^2))
      }
    }
  } else if (d == 2) {
    sigma_mat1 <- matrix(NA, nrow = nrow(r), ncol = nrow(r))
    sigma_mat2 <- matrix(NA, nrow = nrow(r), ncol = nrow(r))
    for (i in seq_len(nrow(r))) {
      for (j in seq_len(nrow(r))) {
        dist1 <- r[i, 1] - r[j, 1]
        dist2 <- r[i, 2] - r[j, 2]
        sigma_mat1[i, j] <- sd^2 * exp(-(dist1^2) / (2 * bandwidth^2))
        sigma_mat2[i, j] <- sd^2 * exp(-(dist2^2) / (2 * bandwidth^2))
      }
    }
    sigma_mat <- sigma_mat1 * sigma_mat2
  }

  return(sigma_mat)
}
