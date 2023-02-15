# highest-density intervals

# true highest-density intervals (from distribution object) -------------------
true_hdi = function(x, mass = 0.95, type = 6) {
  hdi = .hdi_dist(x, mass, type = type)
  matrix(c(hdi$lower, hdi$upper), ncol = 2)
}


# based on hdr.dist_default from {distributional}
# https://github.com/mitchelloharawild/distributional/blob/master/R/default.R
.hdi_dist <- function(dist, mass = 0.95, n = 10000, type = 5, ...) {
  dist_x <- quantile(dist, seq(0.5/n, 1 - 0.5/n, length.out = n))[[1]]
  # Remove duplicate values of dist_x from less continuous distributions
  dist_x <- unique(dist_x)
  dist_y <- density(dist, dist_x)[[1]]
  alpha <- quantile(dist_y, probs = 1 - mass, type = type)


  it <- seq_len(length(dist_y) - 1)
  y_minus_alpha <- dist_y - alpha
  dd <- y_minus_alpha[it + 1] * y_minus_alpha[it]
  index <- it[dd <= 0]
  # unique() removes possible duplicates if sequential dd has same value.
  # More robust approach is required.
  # hdr = unique(vapply(index, FUN.VALUE = numeric(1), function(i) {
  #   uniroot(function(x) density(dist, x)[[1]] - alpha, lower = dist_x[[i]], upper = dist_x[[i + 1]])$root
  # }))
  y0 = y_minus_alpha[index]
  y1 = y_minus_alpha[index + 1]
  x0 = dist_x[index]
  x1 = dist_x[index + 1]
  hdr = unique(x1 - (x1 - x0) / (y1 - y0) * y1)
  # Add boundary values which may exceed the crossing point.
  hdr = c(dist_x[1][dist_y[1] > alpha], hdr, dist_x[length(dist_x)][dist_y[length(dist_y)] > alpha])

  # hdr <- crossing_alpha(alpha, dist_x, dist_y)
  lower_hdr <- seq_along(hdr) %% 2 == 1
  list(lower = hdr[lower_hdr], upper = hdr[!lower_hdr])
}


.hdi_numeric <- function(x, mass = 0.95, n = 4096, density = ggdist::density_bounded, ...) {
  dist_x <- quantile(x, seq(0.5/n, 1 - 0.5/n, length.out = n))
  # Remove duplicate values of dist_x from less continuous distributions
  dist_x <- unique(dist_x)
  d = density(x, n = n)
  dist_y = approx(d$x, d$y, dist_x)$y
  alpha <- quantile(dist_y, probs = 1 - mass)

  it <- seq_len(length(dist_y) - 1)
  y_minus_alpha <- dist_y - alpha
  dd <- y_minus_alpha[it + 1] * y_minus_alpha[it]
  index <- it[dd <= 0]
  # unique() removes possible duplicates if sequential dd has same value.
  y0 = y_minus_alpha[index]
  y1 = y_minus_alpha[index + 1]
  x0 = dist_x[index]
  x1 = dist_x[index + 1]
  hdr = unique(x1 - (x1 - x0) / (y1 - y0) * y1)
  # Add boundary values which may exceed the crossing point.
  hdr = c(dist_x[1][dist_y[1] > alpha], hdr, dist_x[length(dist_x)][dist_y[length(dist_y)] > alpha])

  lower_hdr <- seq_along(hdr) %% 2 == 1
  list(lower = hdr[lower_hdr], upper = hdr[!lower_hdr])
}

# quantile shortest intervals ------------------------------------------------------

# distributional-like implementation of highest-density region
hdru = function(x, mass = .95, options = list()) {
  hdi = .hdi_numeric(x, mass, n = 4096, density = ggdist::density_unbounded)
  matrix(c(hdi$lower, hdi$upper), ncol = 2)
}

# distributional-like implementation of highest-density region, bounded
hdrb = function(x, mass = .95, options = list()) {
  hdi = .hdi_numeric(x, mass, n = 4096, density = ggdist::density_bounded)
  matrix(c(hdi$lower, hdi$upper), ncol = 2)
}

# HDInterval implementation of highest-density interval
hdiu = function(x, mass = .95, options = list()) {
  d = ggdist::density_unbounded(x, n = 4096)
  matrix(HDInterval::hdi(d, credMass = mass, allowSplit = TRUE), ncol = 2)
}

# HDInterval implementation of highest-density interval, bounded
hdib = function(x, mass = .95, options = list()) {
  d = ggdist::density_bounded(x, n = 4096)
  matrix(HDInterval::hdi(d, credMass = mass, allowSplit = TRUE), ncol = 2)
}
