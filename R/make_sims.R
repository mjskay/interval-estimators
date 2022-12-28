# generate sims we will evaluate intervals on

make_sims = function(dist = dist_normal()) {
  expand_grid(
    dist = dist,
    n = c(500, 1000, 10000),
    seed = 1:1000
  ) |>
    mutate(
      sample = mapply(\(dist, seed, n) withr::with_seed(seed, generate(dist, n)), dist, seed, n)
    )
}

