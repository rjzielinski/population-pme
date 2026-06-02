projection_additive_pme <- function(x, f, initial_guess, param_grid, template) {
  require(pme, quietly = TRUE)
  require(pracma, quietly = TRUE)

  n_knots <- nrow(param_grid)
  d <- ncol(param_grid)

  timeval <- x[1]
  x_spatial <- x[-1]
  init_param <- initial_guess[-1]

  time_coefs <- f(timeval)

  coef_mat <- matrix(
    time_coefs,
    nrow = n_knots + d + 1,
    byrow = TRUE
  )

  if (template == "euclidean") {
    coef_s <- t(coef_mat[1:n_knots, ])
    coef_alpha <- t(coef_mat[-(1:n_knots), ])

    embedding_map <- function(t) {
      return_vec <- (coef_s %*% etaFunc(t, param_grid, 4 - d)) +
        (coef_alpha %*% c(1, t))
      return_vec
    }
  } else if (template == "sphere") {
    coef_s <- t(coef_mat[1:n_knots, ])
    coef_alpha <- t(coef_mat[n_knots + 1, , drop = FALSE])

    params_cart <- sph2cart(cbind(param_grid, 1))

    embedding_map <- function(t) {
      t_cart <- pracma::sph2cart(c(t, 1))
      kernel_vals <- sphere_kernel_func(t_cart, params_cart)

      as.vector(
        (coef_s %*% kernel_vals) +
          (coef_alpha %*% matrix(1, ncol = 1))
      )
    }
  }

  obj_fun <- function(t) {
    sum((x_spatial - embedding_map(t))^2)
  }

  nlm_est <- try(
    stats::nlm(
      obj_fun,
      p = init_param
    ),
    silent = TRUE
  )

  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-07)
    nlopt_est <- try(
      nloptr::nloptr(
        x0 = init_param,
        obj_fun,
        opts = opts
      ),
      silent = TRUE
    )
    if (inherits(nlopt_est, "try-error")) {
      return(NULL)
    } else {
      return(c(timeval, nlopt_est$solution))
    }
  } else {
    return(c(timeval, nlm_est$estimate))
  }
}
