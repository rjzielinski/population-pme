scale_spherical <- function(data) {
  theta <- data[, 1]
  phi <- data[, 2]

  # ensure param values are within -pi to pi, or -pi/2 to pi/2
  theta <- sign(theta) * ((abs(theta) / pi) %% 1) * pi
  phi <- sign(phi) * ((abs(phi) / (pi / 2)) %% 1) * (pi / 2)

  cbind(theta, phi)
}
