# equi-tailed intervals

# true equi-tailed intervals (from distribution object) -------------------
true_eti = function(x, mass) {
  ggdist::qi(x, mass)
}


# quantile intervals ------------------------------------------------------

qi = function(x, mass, options = list(type = 7)) {
  lower_prob = (1 - mass)/2
  upper_prob = (1 + mass)/2
  probs = c(lower_prob, upper_prob)

  endpoints = if (isTRUE(options$density)) {
    d = ggdist::density_unbounded(x)
    ggdist::weighted_quantile(d$x, probs, weights = d$y, type = options$type)
  } else {
    quantile(x, probs, type = options$type)
  }
  matrix(endpoints, ncol = 2)
}
