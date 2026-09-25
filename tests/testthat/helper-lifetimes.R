# Helpers for testing the lifetime distribution of simulated ADB trees.
#
# Branches of divided (status 2) and dead (status 0) particles are completed
# lifetimes, but they are NOT a plain sample from the lifetime distribution:
# long lifetimes are more likely to be cut off at the end of the simulation,
# so completed lifetimes are biased towards short values. Alive particles
# (status 1) therefore enter as right-censored observations, and lifetimes are
# tested with a one-sample log-rank test, which accounts for the censoring.

# One row per particle (incl. the root): observed lifetime `time`, `event`
# (TRUE = completed lifetime, FALSE = censored) and `type`.
# censor_alive_at:
#   "tip"           alive particles are censored at their pendant branch length
#                   (valid for the origin simulator: the present is fixed in advance)
#   "last_division" alive particles are censored at the time of the last division
#                   (for the ntaxa simulator, whose pendant branches extend to the
#                   next event and are therefore biased, see note in R/sim_adb_ntaxa.R)
lifetime_data <- function(tree, censor_alive_at = c("tip", "last_division")) {
  censor_alive_at <- match.arg(censor_alive_at)
  ph <- tree@phylo
  n_nodes <- length(ph$tip.label) + ph$Nnode
  root <- length(ph$tip.label) + 1

  edge_length <- numeric(n_nodes)
  edge_length[ph$edge[, 2]] <- ph$edge.length
  edge_length[root] <- ph$root.edge

  data_rows <- match(seq_len(n_nodes), tree@data$node)
  status <- tree@data$status[data_rows]
  type <- tree@data$type[data_rows]

  time <- edge_length
  if (censor_alive_at == "last_division") {
    end <- ape::node.depth.edgelength(ph) + ph$root.edge # time since origin at the end of each lifetime
    birth <- end - edge_length
    alive <- status == 1
    time[alive] <- max(end[status == 2]) - birth[alive]
  }

  data.frame(time = time, event = status != 1, type = type)
}

# One-sample log-rank test of H0: lifetimes of type i ~ Gamma(shape = b[i], scale = a[i]).
#
# Each particle j is observed for a time t_j: its full lifetime if completed
# (event = TRUE, divided or dead), otherwise until it was censored (alive).
# - O = observed number of events = number of completed lifetimes.
# - E = number of events expected under H0 = sum over ALL particles (completed
#   and censored) of the cumulative hazard H(t_j) = -log S(t_j), where S is the
#   survival function of the particle's lifetime distribution. H(t_j) is the
#   expected number of events a particle contributes while observed for t_j.
# Under H0, O - E has mean 0 and variance E, so (O - E)^2 / E ~ chi^2 with 1 df.
# O > E: simulated lifetimes are shorter than H0 (more events than expected);
# O < E: they are longer (e.g. censoring times that are too long).
# This is valid when censoring is independent of the future lifetimes, which is
# why completed lifetimes alone cannot be compared with Gamma(b, a) directly.
logrank_p <- function(lifetimes, a, b) {
  cum_hazard <- -stats::pgamma(lifetimes$time, shape = b[lifetimes$type + 1], scale = a[lifetimes$type + 1],
                               lower.tail = FALSE, log.p = TRUE)
  observed <- sum(lifetimes$event)
  expected <- sum(cum_hazard)
  stats::pchisq((observed - expected)^2 / expected, df = 1, lower.tail = FALSE)
}

# Simulate n trees with sim(), dropping failed simulations (NULL)
simulate_trees <- function(n, sim) {
  trees <- lapply(seq_len(n), function(i) suppressMessages(sim()))
  Filter(Negate(is.null), trees)
}
