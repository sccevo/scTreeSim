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

#' Rate with a time profile (internal)
#'
#' @param base scalar or length-`n_barcodes` base rate
#' @param time a `rate_profile`, e.g. [rate_window()]; `NULL` = constant in time
#' @noRd
barcode_rate <- function(base, time = NULL) {
  if (is.null(time)) time <- rate_constant()
  stopifnot(is.numeric(base), !anyNA(base), all(base >= 0), inherits(time, "rate_profile"))
  structure(list(base = base, time = time), class = "barcode_rate")
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
