library(distributional)
library(dplyr)
library(purrr)

# define the distributions we will simulate data from
dist = list(
  quote(dist_normal(0, 1)),
  quote(dist_normal(0, 2)), # as a check: rmse/sd on N(0,1) and N(0,2) should be the same
  quote(dist_student_t(3, 0, 1)),
  quote(dist_exponential(1)),
  quote(dist_gamma(2, 2)),
  quote(dist_beta(10, 1)),
  quote(dist_beta(10, 3)),
  quote(dist_beta(10, 6)),
  # dist_lognormal(c(0,1,2), c(1,1,0.5)),
  # dist_beta(10, 1),  # as a check: rmse on lower and upper on beta(1,10) should be similar to upper and lower on beta(10,1)
  # dist_beta(c(1,2,3,6,10), 10),
  quote(dist_mixture(dist_normal(), dist_normal(6,2), weights = c(0.6, 0.4)))
  # dist_mixture(dist_normal(), dist_normal(6), weights = c(0.6, 0.4))
)

dist_df = tibble(
  name = map_chr(dist, \(d) gsub("[\\(\\) ,=]+", "_", format(eval(d)))),
  dist = dist
)

sims_df = tibble(
  name = dist_df$name,
  sims = syms(paste0("sims_", dist_df$name))
)

sims_targets = tar_map(
  values = dist_df,
  tar_target(sims, make_sims(dist)),
  names = name
)
