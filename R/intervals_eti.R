# equi-tailed intervals

# true equi-tailed intervals (from distribution object) -------------------
true_eti = function(x, mass) {
  ggdist::qi(x, mass)
}


# quantile intervals ------------------------------------------------------

qi = function(x, mass) {
  ggdist::qi(x, mass)
}
