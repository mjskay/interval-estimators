make_intervals = function(interval_function, sims, mass, type) {
  method = as_label(enquo(interval_function))

  sims |>
    mutate(
      endpoints = lapply(sample, interval_function, mass = mass),
      method = method,
      type = type,
      mass = mass
    ) |>
    select(-sample)
}


# true intervals ----------------------------------------------------------

apply_true_interval_function = function(type, dist, mass) {
  true_interval_function = switch(type,
    eti = true_eti
  )
  true_interval_function(dist, mass)
}

add_true_intervals = function(intervals) {
  intervals |>
    mutate(
      true_endpoints = mapply(apply_true_interval_function, type, dist, mass, SIMPLIFY = FALSE),
    )
}


# coverage ----------------------------------------------------------------

coverage = function(dist, endpoints) {
  sum(mapply(endpoints[,1], endpoints[,2], FUN = \(l, u) cdf(dist, u) - cdf(dist, l)))
}

add_coverage = function(intervals) {
  intervals |>
    mutate(
      coverage = map2_dbl(dist, endpoints, coverage),
    )
}


# percent overlap ---------------------------------------------------------

matrix_to_interval = \(matrix) {
  as.interval(mapply(reals, matrix[,1], matrix[,2], SIMPLIFY = FALSE))
}

percent_overlap = function(x, y) {
  x_interval = matrix_to_interval(x)
  y_interval = matrix_to_interval(y)
  interval_measure(interval_intersection(x, y)) / interval_measure(interval_union(x,y))
}

add_percent_overlap = function(intervals, type) {
  intervals |>
    mutate(
      percent_overlap = map2_dbl(true_endpoints, endpoints, percent_overlap)
    )
}
