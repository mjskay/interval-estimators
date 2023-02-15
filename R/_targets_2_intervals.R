library(rlang)

options_to_string = function(options) {
  gsub("[ =,\\)]+", "_", gsub("list\\(", "", as.character(options)))
}

intervals_df = sims_df |>
  tidyr::expand(
    tidyr::nesting(name, sims),
    tribble(
      ~interval_function, ~type, ~options,

      # equi-tailed interval types
      # quote(qi), "eti", quote(list(type = 4)),
      # quote(qi), "eti", quote(list(type = 5)),
      # quote(qi), "eti", quote(list(type = 6)),
      quote(qi), "eti", quote(list(type = 7)),
      # quote(qi), "eti", quote(list(type = 8)),
      # quote(qi), "eti", quote(list(type = 9)),

      # shortest intervals
      # quote(qsi), "si", quote(list(type = 1)),
      # quote(qsi), "si", quote(list(type = 4)),
      # quote(qsi), "si", quote(list(type = 5)),
      quote(qsi), "si", quote(list(type = 6)),
      # quote(qsi), "si", quote(list(type = 7)),
      # quote(qsi), "si", quote(list(type = 8)),
      # quote(qsi), "si", quote(list(type = 9)),
      # quote(dqsi), "si", quote(list(type = 4)),
      # quote(dqsi), "si", quote(list(type = 5)),
      quote(dqsi), "si", quote(list(type = 5, adjust = 0.5)),
      # quote(dqsi), "si", quote(list(type = 6)),
      # quote(dqsi), "si", quote(list(type = 7)),
      # quote(dqsi), "si", quote(list(type = 8)),
      # quote(dqsi), "si", quote(list(type = 9)),
      quote(spin), "si", quote(list()),
      # highest-density intervals
      quote(hdru), "hdi", quote(list()),
      quote(hdrb), "hdi", quote(list()),
      quote(hdiu), "hdi", quote(list()),
      quote(hdib), "hdi", quote(list())
    ),
    mass = c(.5, .95)
  ) |>
  mutate(
    name = paste0(
      name,
      map_chr(interval_function, as_label), "_",
      options_to_string(options),
      mass * 100
    )
  )

intervals_targets = tar_map(
  values = intervals_df,
  tar_target(intervals_partial_1, make_intervals(interval_function, sims, mass, type, options)),
  tar_target(intervals_partial_2, add_true_intervals(intervals_partial_1)),
  tar_target(intervals_partial_3, add_coverage(intervals_partial_2)),
  tar_target(intervals, add_percent_overlap(intervals_partial_3)),
  names = name
)

intervals_targets_names = map_chr(intervals_targets[[length(intervals_targets)]], \(x) x$settings$name)

# pull out eti intervals --------------------------------------------------

eti_names = paste0("intervals_", intervals_df |> filter(type == "eti") |> pull(name))
eti_intervals_targets = intervals_targets[[length(intervals_targets)]][intervals_targets_names %in% eti_names]

eti_intervals_targets_combined = tar_combine(
  eti_intervals,
  eti_intervals_targets,
  command = dplyr::bind_rows(!!!.x)
)

# pull out shortest intervals --------------------------------------------------

si_names = paste0("intervals_", intervals_df |> filter(type == "si") |> pull(name))
si_intervals_targets = intervals_targets[[length(intervals_targets)]][intervals_targets_names %in% si_names]

si_intervals_targets_combined = tar_combine(
  si_intervals,
  si_intervals_targets,
  command = dplyr::bind_rows(!!!.x)
)

# pull out highest-density intervals --------------------------------------------------

hdi_names = paste0("intervals_", intervals_df |> filter(type == "hdi") |> pull(name))
hdi_intervals_targets = intervals_targets[[length(intervals_targets)]][intervals_targets_names %in% hdi_names]

hdi_intervals_targets_combined = tar_combine(
  hdi_intervals,
  hdi_intervals_targets,
  command = dplyr::bind_rows(!!!.x)
)
