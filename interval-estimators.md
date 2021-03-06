Comparing interval estimators
================
Matthew Kay
2022-07-02

## Setup

Libraries:

``` r
library(distributional)
library(ggdist)
library(bayestestR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(sets)
library(magrittr)

theme_set(theme_ggdist())
dir.create("sims", showWarnings = FALSE)

cache = function(expr, file) {
  if (file.exists(file)) {
    readRDS(file)
  } else {
    saveRDS(expr, file)
    expr
  }
}
```

Some helper functions for pulling out intervals:

``` r
lower = function(x) UseMethod("lower")

lower.hdr = function(x) {
  vapply(vctrs::field(x, "lower"), min, numeric(1))
}

lower.list = function(x) {
  vapply(x, `[[`, numeric(1), "CI_low")
}

upper = function(x) UseMethod("upper")

upper.hdr = function(x) {
  vapply(vctrs::field(x, "upper"), max, numeric(1))
}

upper.list = function(x) {
  vapply(x, `[[`, numeric(1), "CI_high")
}
```

## Distributions to test

The distributions we’ll test on:

``` r
dist = c(
  dist_normal(0, 1:2),  # as a check: rmse/sd on N(0,1) and N(0,2) should be the same
  dist_student_t(3, 0, 1),
  dist_lognormal(c(0,1,2), c(1,1,0.5)),
  dist_beta(10, 1),  # as a check: rmse on lower and upper on beta(1,10) should be similar to upper and lower on beta(10,1)
  dist_beta(c(1,2,3,6,10), 10),
  dist_mixture(dist_normal(), dist_normal(6), weights = c(0.6, 0.4))
)

tibble(dist) |>
  ggplot(aes(xdist = dist)) +
  stat_slab(aes(color = format(dist)), fill = NA, normalize = "panels") +
  facet_wrap(~ family(dist), scales = "free")
```

<img src="interval-estimators_files/figure-gfm/dists-1.png" width="672" />

## Generate data

We’ll generate 100 samples of each sample size from each distribution
using a fixed set of seeds:

``` r
list_of_intervals = \(lower, upper) {
  mapply(reals, lower, upper, SIMPLIFY = FALSE)
}

matrices_to_intervals = \(matrices) {
  lapply(matrices, \(m) as.interval(apply(m, 1, \(x) reals(x[[1]], x[[2]]))))
}

grid = expand.grid(
  dist = dist,
  n = c(500, 1000, 10000),
  #N.B. this is not a sustainable way fo doing this; long-term should not pre-generate the
  #whole grid but do it on-demand and then summarize results from each sample.
  seed = 1:100
) |>
  as_tibble() |>
  mutate(
    sample = mapply(\(dist, seed, n) withr::with_seed(seed, generate(dist, n)), dist, seed, n),
    true_hdci = lapply(dist, ggdist::hdci),
    true_hdci_lower = vapply(true_hdci, min, numeric(1)),
    true_hdci_upper = vapply(true_hdci, max, numeric(1)),
    true_hdci = list_of_intervals(true_hdci_lower, true_hdci_upper),
    
    true_hdi = distributional::hdr(dist),
    true_hdi = matrices_to_intervals(mapply(cbind, 
      vctrs::field(true_hdi, "lower"), vctrs::field(true_hdi, "upper"), SIMPLIFY = FALSE
    )),
    # for non-mixture distributions the true hdi is the hdci
    true_hdi = ifelse(family(dist) == "mixture", true_hdi, true_hdci),
    true_hdi_lower = vapply(true_hdi, min, numeric(1)),
    true_hdi_upper = vapply(true_hdi, max, numeric(1))
  ) |>
  cache("sims/grid.rds")
```

## Calculate intervals

Calculate HDCIs

``` r
hdci_df = grid |>
  mutate(
    est_interval = lapply(sample, bayestestR::hdi),
    est_lower = lower(est_interval),
    est_upper = upper(est_interval),
    method = "hdci",
    type = "hdci",
    est_interval = list_of_intervals(est_lower, est_upper)
  ) |>
  rename(
    true_interval = true_hdci,
    lower = true_hdci_lower,
    upper = true_hdci_upper
  ) |>
  select(-sample, -true_hdi, -true_hdi_lower, -true_hdi_upper) |>
  cache("sims/hdci_df.rds")
```

Calculate HDCIs using the implementation as in HDInterval (via ggdist):

``` r
hdci2_df = grid |>
  mutate(
    est_interval = lapply(sample, ggdist::hdci),
    est_lower = vapply(est_interval, min, numeric(1)),
    est_upper = vapply(est_interval, max, numeric(1)),
    method = "hdci2",
    type = "hdci",    
    est_interval = list_of_intervals(est_lower, est_upper)
  ) |>
  rename(
    true_interval = true_hdci,
    lower = true_hdci_lower,
    upper = true_hdci_upper
  ) |>
  select(-sample, -true_hdi, -true_hdi_lower, -true_hdi_upper) |>
  cache("sims/hdci2_df.rds")
```

Calculate HDIs. Note we have to take the lower limit of the lower
interval (resp. upper limit of upper interval) as there may be multiple
intervals per sample:

``` r
hdi_df = grid |>
  mutate(
    est_interval = lapply(sample, ggdist::hdi),
    est_lower = vapply(est_interval, min, numeric(1)),
    est_upper = vapply(est_interval, max, numeric(1)),
    method = "hdi",
    type = "hdi",
    est_interval = matrices_to_intervals(est_interval)
  ) |>
  rename(
    true_interval = true_hdi,
    lower = true_hdi_lower,
    upper = true_hdi_upper
  ) |>
  select(-sample, -true_hdci, -true_hdci_lower, -true_hdci_upper) |>
  cache("sims/hdi_df.rds")
```

Calculate SPIs:

``` r
spi_df = grid |>
  mutate(
    est_interval = lapply(sample, bayestestR::spi),
    est_lower = lower(est_interval),
    est_upper = upper(est_interval),
    method = "spi",
    type = "hdci",    
    est_interval = list_of_intervals(est_lower, est_upper)    
  ) |>
  rename(
    true_interval = true_hdci,
    lower = true_hdci_lower,
    upper = true_hdci_upper
  ) |>
  select(-sample, -true_hdi, -true_hdi_lower, -true_hdi_upper) |>
  cache("sims/spi_df.rds")
```

Calculate QIs:

``` r
qi_df = grid |>
  mutate(
    est_interval = lapply(sample, bayestestR::eti),
    est_lower = lower(est_interval),
    est_upper = upper(est_interval),
    method = "qi",
    type = "qi",    
    lower = quantile(dist, 0.025),
    upper = quantile(dist, 0.975),
    true_interval = list_of_intervals(lower, upper),
    est_interval = list_of_intervals(est_lower, est_upper)
  ) |>
  select(-sample, -true_hdi, -true_hdi_lower, -true_hdi_upper, -true_hdci, -true_hdci_lower, -true_hdci_upper) |>
  cache("sims/qi_df.rds")
```

Calculate HDRs using method from hdrcde:

``` r
hdrcde_df = grid |>
  mutate(
    est_interval = lapply(sample, \(x) matrix(hdrcde::hdr(x, prob = 95)$hdr, byrow = TRUE, ncol = 2)),
    est_lower = vapply(est_interval, min, numeric(1)),
    est_upper = vapply(est_interval, max, numeric(1)),
    method = "hdrcde",
    type = "hdi",
    est_interval = matrices_to_intervals(est_interval)
  ) |>
  rename(
    true_interval = true_hdi,
    lower = true_hdi_lower,
    upper = true_hdi_upper
  ) |>
  select(-sample, -true_hdci, -true_hdci_lower, -true_hdci_upper) |>
  cache("sims/hdrcde_df.rds")
```

Calculate HDRs using method from distributional:

``` r
hdr_df = grid |>
  mutate(
    est_interval = distributional::hdr(dist_sample(sample)),
    est_lower = lower(est_interval),
    est_upper = upper(est_interval),
    method = "hdr",
    type = "hdi",    
    est_interval = matrices_to_intervals(mapply(cbind, 
      vctrs::field(est_interval, "lower"), vctrs::field(est_interval, "upper"), SIMPLIFY = FALSE
    ))
  ) |>
  rename(
    true_interval = true_hdi,
    lower = true_hdi_lower,
    upper = true_hdi_upper
  ) |>
  select(-sample, -true_hdci, -true_hdci_lower, -true_hdci_upper) |>
  cache("sims/hdr_df.rds")
```

Shortest interval using quantile estimator:

``` r
quantile_shortest_interval = \(x, mass = 0.95, ...) {
  interval_width = \(p) quantile(x, p + mass, ...) - quantile(x, p, ...)
  p = optimize(interval_width, c(0, 1 - mass))$minimum
  as.vector(quantile(x, c(p, p + mass), ...))
}

qsi_df = grid |>
  mutate(
    est_interval = lapply(sample, quantile_shortest_interval, type = 8),
    est_lower = vapply(est_interval, min, numeric(1)),
    est_upper = vapply(est_interval, max, numeric(1)),
    method = "qsi",
    type = "hdci",
    est_interval = list_of_intervals(est_lower, est_upper)
  ) |>
  rename(
    true_interval = true_hdci,
    lower = true_hdci_lower,
    upper = true_hdci_upper
  ) |>
  select(-sample, -true_hdi, -true_hdi_lower, -true_hdi_upper) |>
  cache("sims/qsi_df.rds")
```

## Comparing methods

As a simple comparison, we’ll look at the root-mean-squared error of
both interval limits compared to their true values. We scale the RMSE by
the standard deviation of the distribution just for comparison purposes
(as some of these distributions have quite different scales):

``` r
dist_names = format(dist)

all_df = bind_rows(
    hdi_df,
    hdci2_df,
    hdci_df,
    hdr_df,
    hdrcde_df,
    qi_df,
    spi_df,
    qsi_df
  ) |>
  # enforce a sensible display order on dist names
  mutate(dist_name = factor(format(dist), levels = rev(dist_names)))

scale_shape_method = scale_shape_manual(values = c(scales::shape_pal()(6), c(1, 2)))

all_df |>
  group_by(method, dist, dist_name, n) |>
  summarise(
    lower = sqrt(mean((est_lower - lower)^2)),
    upper = sqrt(mean((est_upper - upper)^2)), 
    .groups = "drop"
  ) |>
  pivot_longer(c(lower, upper), names_to = "limit", values_to = "rmse") |>
  ggplot(aes(x = rmse / sqrt(variance(dist)), y = dist_name, color = method, shape = method)) +
  geom_point(alpha = 0.75, size = 3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, alpha = 0.25) +
  facet_grid(limit ~ n) +
  scale_shape_method
```

<img src="interval-estimators_files/figure-gfm/limit_plot-1.png" width="672" />

The “true” instance of an interval (actually set of intervals) for a
given distribution and interval type (e.g. quantile interval,
highest-density interval (set), or shortest interval) is some subset of
![\\mathbb{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D "\mathbb{R}"),
![I](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I "I").
Each interval estimator yields an estimate
![\\hat I](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%20I "\hat I"),
expressed as a set of disjoint intervals
![\\left\\{\[L_1, U_1\], \\dots\\right\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%5C%7B%5BL_1%2C%20U_1%5D%2C%20%5Cdots%5Cright%5C%7D "\left\{[L_1, U_1], \dots\right\}").
Given the sum of the widths of those intervals (i.e. the Lesbegue
measure of the set of intervals,
![l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;l "l")),
the percent overlap of the estimated interval set and the true interval
set is:

![
\\textrm{percent overlap}(I, \\hat{I}) = \\frac{l(I \\cap \\hat{I})}{l(I \\cup \\hat{I})}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Ctextrm%7Bpercent%20overlap%7D%28I%2C%20%5Chat%7BI%7D%29%20%3D%20%5Cfrac%7Bl%28I%20%5Ccap%20%5Chat%7BI%7D%29%7D%7Bl%28I%20%5Ccup%20%5Chat%7BI%7D%29%7D%0A "
\textrm{percent overlap}(I, \hat{I}) = \frac{l(I \cap \hat{I})}{l(I \cup \hat{I})}
")

This has the advantage (compared to just comparing lower and upper
endpoints) that it generalizes to estimators that return multiple
intervals.

We can calculate the percent overlap for all interval estimates with
their true interval:

``` r
percent_overlap = \(x, y) {
  map2_dbl(x, y, \(x, y) interval_measure(interval_intersection(x, y)) / interval_measure(interval_union(x,y)))
}

overlap_df = all_df |>
  mutate(percent_overlap = percent_overlap(true_interval, est_interval)) |>
  cache("sims/overlap_df.rds")
```

And compare:

``` r
overlap_df |>
  group_by(method, dist, dist_name, n) |>
  summarise(
    mean_percent_overlap = mean(percent_overlap),
    .groups = "drop"
  ) |>
  ggplot(aes(x = mean_percent_overlap, y = dist_name, color = method, shape = method)) +
  geom_point(alpha = 0.75, size = 3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, alpha = 0.25) +
  facet_grid( ~ n) +
  scale_shape_method
```

<img src="interval-estimators_files/figure-gfm/percent_overlap_plot-1.png" width="672" />

And finally, coverage:

``` r
coverage = \(dist, interval) {
  sum(vapply(unclass(interval), \(i) cdf(dist, i$r) - cdf(dist, i$l), numeric(1)))
}

coverage_df = all_df |>
  mutate(
    coverage = map2_dbl(dist, est_interval, coverage)
  ) |>
  cache("sims/coverage_df.rds")
```

``` r
coverage_df |>
  group_by(method, dist, dist_name, n) |>
  summarise(
    mean_coverage = mean(coverage),
    .groups = "drop"
  ) |>
  ggplot(aes(x = mean_coverage, y = dist_name, color = method, shape = method)) +
  geom_point(alpha = 0.75, size = 3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0.95, alpha = 0.25) +
  facet_grid(n ~ .) +
  scale_shape_method 
```

<img src="interval-estimators_files/figure-gfm/coverage_plot-1.png" width="768" />

## Variance

``` r
all_df |> 
  group_by(method, type, dist, dist_name, n) |> 
  summarise(sd = sd(est_lower)) |>
  ggplot(aes(x = sd/sqrt(variance(dist)), y = dist_name, color = method, shape = method)) +
  geom_point(alpha = 0.75, size = 3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, alpha = 0.25) +
  facet_grid(n ~ .) +
  scale_shape_method 
```

    ## `summarise()` has grouped output by 'method', 'type', 'dist', 'dist_name'. You
    ## can override using the `.groups` argument.

<img src="interval-estimators_files/figure-gfm/variance_plot-1.png" width="768" />

Also worth looking at is the variance of the interval endpoints. I
haven’t gone far with that yet except to show a comparison of QIs to
some HDCI implementations:

``` r
set.seed(1234)
df = data.frame(
  x = 1:100,
  y = rnorm(1000000, 0)
) 

bind_rows(
  df |>
    group_by(x) |>
    median_qi(y, .width = .8) %T>%
    {cat("qi:\t", sd(.$.lower), "\n")} %>%
    mutate(method = "qi"),
  
  df |>
    group_by(x) |>
    median_hdci(y, .width = .8) %T>%
    {cat("hdci:\t", sd(.$.lower), "\n")} %>%
    mutate(method = "hdci"),

  df |>
    group_by(x) |>
    summarise(
      qsi = list(quantile_shortest_interval(y, type = 8, mass = 0.8)),
      .lower = qsi[[1]][[1]],
      .upper = qsi[[1]][[2]],
      y = median(y),
      method = "qsi"
    ) %T>%
    {cat("qsi:\t", sd(.$.lower), "\n")},
  
  df |>
    group_by(x) |>
    summarise(
      spi = list(bayestestR::spi(y, ci = .8)),
      .lower = spi[[1]]$CI_low,
      .upper = spi[[1]]$CI_high,
      y = median(y),
      method = "spi"
    ) %T>%
    {cat("spi:\t", sd(.$.lower), "\n")},

  df |>
    group_by(x) |>
    summarise(
      hdi = list(distributional::hdr(dist_sample(list(y)), size = 80)),
      .lower = lower(hdi[[1]]),
      .upper = upper(hdi[[1]]),
      y = median(y),
      method = "hdr"
    ) %T>%
    {cat("hdr:\t", sd(.$.lower), "\n")}
) |>
  ggplot(aes(x = x, y = y, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(fill = "gray65") +
  ggtitle("qi") +
  ylim(-4, 4) +
  facet_wrap(~ method)
```

    ## qi:   0.01716029 
    ## hdci:     0.04776004 
    ## qsi:  0.0418678 
    ## spi:  0.04773497 
    ## hdr:  0.0214488

<img src="interval-estimators_files/figure-gfm/variance_ribbons-1.png" width="672" />
