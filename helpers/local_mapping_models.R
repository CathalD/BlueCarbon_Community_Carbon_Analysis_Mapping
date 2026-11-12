# Distance models used by both local and Janousek-prior mapping.

distance_matrix <- function(x, y) {
  sqrt(outer(x[, 1], y[, 1], "-")^2 + outer(x[, 2], y[, 2], "-")^2)
}

# Matérn distance scale: at smoothness 0.5, correlation = exp(-distance / range_m).
matern_correlation <- function(distances, range_m, smoothness) {
  scaled <- distances / range_m
  positive <- scaled > 0
  z <- scaled[positive]
  correlation <- scaled
  correlation[] <- 1
  correlation[positive] <- pmin(1, exp(
    (1 - smoothness) * log(2) - lgamma(smoothness) + smoothness * log(z) +
      log(besselK(z, smoothness, expon.scaled = TRUE)) - z
  ))
  correlation
}

idw_predict <- function(distances, log_stock, power, neighbours) {
  vapply(seq_len(ncol(distances)), function(j) {
    d <- distances[, j]
    coincident <- d == 0
    if (any(coincident)) return(mean(log_stock[coincident]))
    nearest <- head(order(d), min(neighbours, length(d)))
    # Scaling all distances by their minimum leaves normalized weights unchanged.
    weights <- (min(d[nearest]) / d[nearest])^power
    weighted.mean(log_stock[nearest], weights)
  }, numeric(1))
}

matern_predict <- function(core_distances, target_distances, log_stock,
                           range_m, smoothness, nugget_fraction) {
  # Work with unit total variance first. The nugget contributes to measured
  # core variance, but not to the latent sediment stock surface.
  process_fraction <- 1 - nugget_fraction
  observation_cov <- process_fraction * matern_correlation(core_distances, range_m, smoothness)
  diag(observation_cov) <- diag(observation_cov) + nugget_fraction
  lower <- t(chol(observation_cov))
  intercept <- as.vector(forwardsolve(lower, rep(1, length(log_stock))))
  whitened_stock <- as.vector(forwardsolve(lower, log_stock))
  intercept_precision <- sum(intercept^2)
  log_mean <- sum(intercept * whitened_stock) / intercept_precision
  residual <- whitened_stock - intercept * log_mean
  # Residual variance estimate at the selected range and smoothness.
  sill <- sum(residual^2) / (length(log_stock) - 1)
  if (all(log_stock == log_stock[1])) sill <- 0

  cross_cov <- process_fraction * matern_correlation(target_distances, range_m, smoothness)
  conditioned <- forwardsolve(lower, cross_cov)
  mean_term <- (1 - drop(crossprod(intercept, conditioned))) / sqrt(intercept_precision)
  list(
    log_mean = log_mean + drop(crossprod(residual, conditioned)),
    log_variance = sill * pmax(0, process_fraction - colSums(conditioned^2) + mean_term^2),
    process_log_variance = sill * process_fraction,
    observation_log_sd = sqrt(sill * nugget_fraction),
    conditioned = sqrt(sill) * conditioned,
    mean_term = sqrt(sill) * mean_term
  )
}

log_stock_moments <- function(log_mean, log_variance) {
  mean_stock <- exp(log_mean + log_variance / 2)
  data.frame(mean_kg_c_m2 = mean_stock, sd_kg_c_m2 = mean_stock * sqrt(expm1(log_variance)))
}

lognormal_total_sd <- function(mean_stock, area_m2, xy, range_m, smoothness,
                               process_log_variance, conditioned, mean_term = NULL) {
  # A joint total needs every pair of cells. Calculate in row blocks and
  # leave the total SD unavailable above 10,000 cells (100 million pairs).
  n <- length(mean_stock)
  if (n > 10000L) return(NA_real_)
  weights_mg <- mean_stock * area_m2 / 1000
  total_variance <- 0
  for (start in seq.int(1L, n, by = 256L)) {
    rows <- start:min(n, start + 255L)
    block_distances <- distance_matrix(xy[rows, , drop = FALSE], xy)
    log_covariance <- process_log_variance *
      matern_correlation(block_distances, range_m, smoothness) -
      crossprod(conditioned[, rows, drop = FALSE], conditioned)
    if (!is.null(mean_term)) log_covariance <- log_covariance + outer(mean_term[rows], mean_term)
    total_variance <- total_variance + sum(outer(weights_mg[rows], weights_mg) * expm1(log_covariance))
  }
  sqrt(max(0, total_variance))
}
