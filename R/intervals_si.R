# shortest intervals

# true shortest intervals (from distribution object) -------------------
true_si = function(x, mass) {
  ggdist::hdci(x, mass)
}


# quantile shortest intervals ------------------------------------------------------

qsi = function(x, mass = .95, options = list(type = 6)) {
  type = options$type
  quantile_fun = ggdist::weighted_quantile_fun(x, type = type)
  p = optimise(function(p) quantile_fun(p + mass) - quantile_fun(p), lower = 0, upper = 1 - mass)$minimum
  endpoints = quantile_fun(c(p, p + mass))
  matrix(endpoints, ncol = 2)
}

dqsi = function(x, mass = .95, options = list(type = 6, adjust = 1), density = ggdist::density_bounded) {
  type = options$type
  adjust = options$adjust %||% 1
  d = density(x, bandwidth = "SJ", adjust = adjust)
  quantile_fun = ggdist::weighted_quantile_fun(d$x, weights = d$y, type = type)
  p = optimise(function(p) quantile_fun(p + mass) - quantile_fun(p), lower = 0, upper = 1 - mass)$minimum
  endpoints = quantile_fun(c(p, p + mass))
  matrix(endpoints, ncol = 2)
}

spin = function(x, mass = .95, options = list()) {
  endpoints = bayestestR::spi(x, mass)
  matrix(c(endpoints$CI_low, endpoints$CI_high), ncol = 2)
}
