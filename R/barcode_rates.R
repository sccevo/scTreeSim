# Time profiles for the rates of the generic barcode simulator.
# A profile is a multiplier g(t) >= 0 of a base rate, as a function of the time since the origin.
# It carries its breakpoints (the simulator splits branches there) and
# max(a, b), an upper bound of g on a piece [a, b] containing no breakpoint
# (needed for thinning).
# All profiles except rate_fn() also carry `pexp`, the closed form
# g(t) = values[findInterval(t, breaks) + 1] * exp(-decay * t) (piecewise constant times an exponential
# decay), which the simulator integrates over a branch instead of thinning.

.rate_profile <- function(g, breaks, max, pexp = NULL) {
  structure(list(g = g, breaks = breaks, max = max, pexp = pexp), class = "rate_profile")
}

.pexp <- function(breaks = numeric(0), values = 1, decay = 0) {
  list(breaks = breaks, values = values, decay = decay)
}

#' Time profiles of rates
#'
#' Multipliers `g(t) >= 0` of a base rate (see [barcode_rate()]) as a function of the
#' time `t` since the origin of the tree, for the generic barcode simulator
#' [sim_barcode_generic()]. Only the integral of the rate over a branch matters for the state at
#' its end, which the simulator computes in closed form for all profiles except `rate_fn()` (also inside
#' `rate_product()`). Those are simulated by thinning instead, which is slower: branches are cut at the
#' breakpoints of a profile, and candidate events are accepted with probability `g(t)` divided by its
#' maximum on the piece.
#'
#' * `rate_constant()`: multiplier 1 at all times.
#' * `rate_window(start, end)`: multiplier 1 on `[start, end)`, 0 elsewhere.
#' * `rate_exp_decay(half_life)`: multiplier `2^(-t / half_life)`, i.e. 1 at the origin and halved
#'   every `half_life` time units.
#' * `rate_piecewise(breaks, values)`: multiplier `values[k]` between `breaks[k - 1]` and
#'   `breaks[k]`; `values` has one more entry than `breaks`.
#' * `rate_fn(f, max)`: an arbitrary function `f(t) >= 0` with a global upper bound `max`. The bound is
#'   checked at the simulated candidate times; a value above it is an error, since it would bias the rate.
#' * `rate_product(...)`: the product of profiles, e.g. decay inside an editing window,
#'   `rate_product(rate_window(1, 2), rate_exp_decay(0.5))`. The decay is measured from the origin.
#'
#' @param start,end start and end of the window (time since the origin), `start <= end`
#' @param half_life time after which the multiplier is halved, positive
#' @param breaks strictly increasing breakpoints (time since the origin)
#' @param values non-negative multipliers, one per interval defined by `breaks`
#' @param f function of a single time returning a non-negative number
#' @param max upper bound of `f`
#' @param ... profiles created by the functions above
#'
#' @return an object of class `rate_profile`, to be used as the `time` argument of [barcode_rate()]
#' @family barcode rates
#' @name rate_profiles
NULL

#' @rdname rate_profiles
#' @export
rate_constant <- function() {
  .rate_profile(function(t) 1, numeric(0), function(a, b) 1, .pexp())
}

#' @rdname rate_profiles
#' @export
rate_window <- function(start, end) {
  stopifnot(is.numeric(start), is.numeric(end), length(start) == 1, length(end) == 1, start <= end)
  g <- function(t) as.numeric(t >= start & t < end)
  .rate_profile(g, c(start, end), function(a, b) g((a + b) / 2), .pexp(c(start, end), c(0, 1, 0)))
}

#' @rdname rate_profiles
#' @export
rate_exp_decay <- function(half_life) {
  stopifnot(is.numeric(half_life), length(half_life) == 1, half_life > 0)
  g <- function(t) 2^(-t / half_life)
  .rate_profile(g, numeric(0), function(a, b) g(a), .pexp(decay = log(2) / half_life))
}

#' @rdname rate_profiles
#' @export
rate_piecewise <- function(breaks, values) {
  stopifnot(
    is.numeric(breaks), !is.unsorted(breaks, strictly = TRUE),
    is.numeric(values), length(values) == length(breaks) + 1, all(values >= 0)
  )
  g <- function(t) values[findInterval(t, breaks) + 1L]
  .rate_profile(g, breaks, function(a, b) g((a + b) / 2), .pexp(breaks, values))
}

#' @rdname rate_profiles
#' @export
rate_fn <- function(f, max) {
  stopifnot(is.function(f), is.numeric(max), length(max) == 1, max >= 0)
  .rate_profile(f, numeric(0), function(a, b) max)
}

#' @rdname rate_profiles
#' @export
rate_product <- function(...) {
  profiles <- list(...)
  stopifnot(length(profiles) >= 1, all(vapply(profiles, inherits, logical(1), "rate_profile")))
  # breakpoints are pooled, so on a piece between them every factor is within its own bound
  breaks <- sort(unique(unlist(lapply(profiles, `[[`, "breaks"))))
  # closed form if every factor has one: the product of the factors on each interval
  # [breaks[k - 1], breaks[k]), evaluated at its left end, and the sum of the decays
  specs <- lapply(profiles, `[[`, "pexp")
  pexp <- if (!any(vapply(specs, is.null, logical(1)))) {
    left <- c(-Inf, sort(unique(unlist(lapply(specs, `[[`, "breaks")))))
    .pexp(left[-1],
          Reduce(`*`, lapply(specs, function(s) s$values[findInterval(left, s$breaks) + 1L])),
          sum(vapply(specs, `[[`, numeric(1), "decay")))
  }
  .rate_profile(
    function(t) Reduce(`*`, lapply(profiles, function(p) p$g(t))),
    breaks,
    function(a, b) prod(vapply(profiles, function(p) p$max(a, b), numeric(1))),
    pexp
  )
}

#' Rate with a time profile and type multipliers
#'
#' Describes an editing or silencing rate (or a dropout probability) of the generic barcode simulator
#' [sim_barcode_generic()] that varies in time and/or between cell types. The rate on a branch of type
#' `k` at time `t` since the origin is `base * by_type[k + 1] * time$g(t)`. The type of a branch is the
#' `type` of the node it ends in (the root edge takes the root's type), which is exact for trees
#' that keep their type changes (`collapse = FALSE`) and approximate otherwise.
#'
#' A dropout probability is evaluated at the tips only, so it takes `by_type` but no `time` profile,
#' and `base * by_type` must not exceed 1 for the types at the tips.
#'
#' @param base scalar or length-`n_barcodes` base rate
#' @param time a `rate_profile`, e.g. [rate_window()]; `NULL` (default) = constant in time
#' @param by_type non-negative multipliers indexed by type + 1 (one entry per type 0, 1, ...);
#'   `NULL` (default) = same for all types. Entries may be `NA` for types that are never used,
#'   e.g. types not occurring at tips for dropout.
#'
#' @return an object of class `barcode_rate`, for the arguments `edit_rate`, `silencing_rate` and
#'   `dropout_prob` of [sim_barcode_generic()]
#' @family barcode rates
#' @export
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

# Integral of a closed-form profile (`pexp`) over each branch [t0, t1] (vectorised over branches)
.profile_integral <- function(pexp, t0, t1) {
  lo <- c(-Inf, pexp$breaks)
  hi <- c(pexp$breaks, Inf)
  total <- numeric(length(t0))
  for (k in seq_along(pexp$values)) {
    if (pexp$values[k] == 0) next
    a <- pmax(t0, lo[k])
    b <- pmin(t1, hi[k])
    inside <- b > a
    a <- a[inside]
    b <- b[inside]
    piece <- if (pexp$decay == 0) b - a else exp(-pexp$decay * a) * -expm1(-pexp$decay * (b - a)) / pexp$decay
    total[inside] <- total[inside] + pexp$values[k] * piece
  }
  total
}

# A rate (editing or silencing) on every branch, as passed to sim_barcodes_cpp(): the branch ending in
# node i spans [t0[i], t1[i]] (time since the origin), and its type multiplies the rate by mult[i].
# Closed-form profiles give the integrated rate of each branch; otherwise each branch is cut at the
# breakpoints of the profile into pieces with their bounds, for thinning.
.rate_process <- function(rate, t0, t1, mult) {
  profile <- rate$time
  if (!is.null(profile$pexp)) {
    return(list(scale = mult * .profile_integral(profile$pexp, t0, t1), g = NULL))
  }
  pieces <- lapply(seq_along(t0), function(i) {
    cuts <- sort(unique(c(t0[i], t1[i], profile$breaks)))
    cuts <- cuts[cuts >= t0[i] & cuts <= t1[i]]
    a <- cuts[-length(cuts)]
    b <- cuts[-1]
    list(start = a, end = b, max = as.numeric(mapply(profile$max, a, b)))
  })
  n_pieces <- vapply(pieces, function(p) length(p$start), integer(1))
  list(scale = mult, g = profile$g, ptr = c(0L, cumsum(n_pieces)),
       start = unlist(lapply(pieces, `[[`, "start")), end = unlist(lapply(pieces, `[[`, "end")),
       max = unlist(lapply(pieces, `[[`, "max")))
}
