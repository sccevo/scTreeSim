# Time profiles for the rates of the generic barcode simulator (internal).
# A profile is a multiplier g(t) >= 0 of a base rate, as a function of the time since the origin.
# It carries its breakpoints (the simulator splits branches there) and
# max(a, b), an upper bound of g on a piece [a, b] containing no breakpoint
# (needed for thinning).

.rate_profile <- function(g, breaks, max) {
  structure(list(g = g, breaks = breaks, max = max), class = "rate_profile")
}

# multiplier 1 at all times
rate_constant <- function() {
  .rate_profile(function(t) 1, numeric(0), function(a, b) 1)
}

# multiplier 1 on [start, end), 0 elsewhere (time since origin)
rate_window <- function(start, end) {
  stopifnot(is.numeric(start), is.numeric(end), length(start) == 1, length(end) == 1, start <= end)
  g <- function(t) as.numeric(t >= start & t < end)
  .rate_profile(g, c(start, end), function(a, b) g((a + b) / 2))
}

# multiplier 2^(-t / half_life), i.e. exponentially decaying activity
rate_exp_decay <- function(half_life) {
  stopifnot(is.numeric(half_life), length(half_life) == 1, half_life > 0)
  g <- function(t) 2^(-t / half_life)
  .rate_profile(g, numeric(0), function(a, b) g(a))
}

# multiplier values[k] between breaks[k - 1] and breaks[k]; length(values) = length(breaks) + 1
rate_piecewise <- function(breaks, values) {
  stopifnot(
    is.numeric(breaks), !is.unsorted(breaks, strictly = TRUE),
    is.numeric(values), length(values) == length(breaks) + 1, all(values >= 0)
  )
  g <- function(t) values[findInterval(t, breaks) + 1L]
  .rate_profile(g, breaks, function(a, b) g((a + b) / 2))
}

# arbitrary multiplier f(t) >= 0 with a global upper bound `max`;
# values above the bound are an error (they would bias the thinning)
rate_fn <- function(f, max) {
  stopifnot(is.function(f), is.numeric(max), length(max) == 1, max >= 0)
  .rate_profile(f, numeric(0), function(a, b) max)
}

# product of profiles, e.g. exponential decay inside a window:
# rate_product(rate_window(1, 2), rate_exp_decay(0.5)).
# Breakpoints are pooled, so on a piece between them every factor is within its own bound
rate_product <- function(...) {
  profiles <- list(...)
  stopifnot(length(profiles) >= 1, all(vapply(profiles, inherits, logical(1), "rate_profile")))
  .rate_profile(
    function(t) Reduce(`*`, lapply(profiles, function(p) p$g(t))),
    sort(unique(unlist(lapply(profiles, `[[`, "breaks")))),
    function(a, b) prod(vapply(profiles, function(p) p$max(a, b), numeric(1)))
  )
}

#' Rate with a time profile and type multipliers (internal)
#'
#' The rate on a branch of type `k` at time `t` is `base x by_type[k + 1] x time$g(t)`.
#'
#' @param base scalar or length-`n_barcodes` base rate
#' @param time a `rate_profile`, e.g. [rate_window()]; `NULL` = constant in time
#' @param by_type non-negative multipliers indexed by type + 1; `NULL` = same for all types.
#'   Entries may be `NA` for types that are never used (e.g. internal-only types for dropout)
#' @noRd
barcode_rate <- function(base, time = NULL, by_type = NULL) {
  time_given <- !is.null(time)
  if (!time_given) time <- rate_constant()
  stopifnot(
    is.numeric(base), !anyNA(base), all(base >= 0), inherits(time, "rate_profile"),
    is.null(by_type) || (is.numeric(by_type) && length(by_type) >= 1 && all(by_type >= 0, na.rm = TRUE))
  )
  structure(list(base = base, time = time, time_given = time_given, by_type = by_type), class = "barcode_rate")
}

# normalise a plain number / vector or a barcode_rate to a barcode_rate
.as_barcode_rate <- function(x, n_barcodes, name) {
  if (!inherits(x, "barcode_rate")) x <- barcode_rate(x)
  if (!length(x$base) %in% c(1, n_barcodes)) {
    stop("`", name, "` must have length 1 or `n_barcodes`.", call. = FALSE)
  }
  if (length(x$base) == 1) x$base <- rep(x$base, n_barcodes)
  x
}

# types of the tree nodes, checked against the `by_type` multipliers of the rates
# (NULL if none of the rates depends on the type). Edit and silencing act on every branch,
# dropout only at the tips, so only the tip types need a (non-NA) dropout multiplier.
.check_node_types <- function(tree_df, edit_rate, silencing_rate, dropout_prob, tip_nodes) {
  rates <- list(edit_rate = edit_rate, silencing_rate = silencing_rate, dropout_prob = dropout_prob)
  if (all(vapply(rates, function(r) is.null(r$by_type), logical(1)))) return(NULL)
  if (!"type" %in% names(tree_df) || anyNA(tree_df$type) || any(tree_df$type < 0) ||
      any(tree_df$type != round(tree_df$type))) {
    stop("`by_type` needs a `type` (0, 1, ...) for every node of the tree.", call. = FALSE)
  }
  node_type <- stats::setNames(tree_df$type, tree_df$node)
  used_types <- list(
    edit_rate = unique(node_type), silencing_rate = unique(node_type),
    dropout_prob = unique(node_type[as.character(tip_nodes)])
  )
  for (name in names(rates)) {
    by_type <- rates[[name]]$by_type
    if (is.null(by_type)) next
    if (length(by_type) <= max(node_type)) {
      stop("`by_type` of `", name, "` must have an entry for every type 0, ..., ", max(node_type), ".", call. = FALSE)
    }
    if (anyNA(by_type[used_types[[name]] + 1])) {
      stop("`by_type` of `", name, "` is NA for a type that occurs where it is used",
           if (name == "dropout_prob") " (the tips)." else " (any node).", call. = FALSE)
    }
  }
  node_type
}

# multiplier of a rate on a branch ending in a node of the given type
.type_multiplier <- function(rate, type) {
  if (is.null(rate$by_type)) 1 else rate$by_type[type + 1]
}
