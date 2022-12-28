library(rlang)

intervals_df = sims_df |>
  tidyr::expand(
    tidyr::nesting(name, sims),
    tribble(
      ~interval_function, ~type,
      quote(qi), "eti"
    ),
    mass = c(.5, .95)
  ) |>
  mutate(
    name = paste0(name, map_chr(interval_function, as_label), "_", mass * 100)
  )

intervals_targets = tar_map(
  values = intervals_df,
  tar_target(intervals_partial_1, make_intervals(interval_function, sims, mass, type)),
  tar_target(intervals_partial_2, add_true_intervals(intervals_partial_1)),
  tar_target(intervals_partial_3, add_coverage(intervals_partial_2)),
  tar_target(intervals_partial_4, add_percent_overlap(intervals_partial_3)),
  tar_target(intervals, add_true_intervals(intervals_partial_4)),
  names = name
)

intervals_targets_combined = tar_combine(
  intervals,
  intervals_targets[[length(intervals_targets)]],
  command = dplyr::bind_rows(!!!.x)
)
