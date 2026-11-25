#' Prior distribution of a parameter
#'
#' Describes the distribution that a simulated parameter is drawn from, optionally
#' truncated to an interval, e.g. the parameters of the genes in
#' \code{\link{sim_gene_params}} and the dispersion of the counts in
#' \code{\link{sim_count_params}}. \code{\link{rsim_prior}} draws from it.
#'
#' The parameters of each family are
#' \itemize{
#'   \item \code{"exponential"}: \code{rate}
#'   \item \code{"lognormal"}: \code{meanlog}, \code{sdlog}
#'   \item \code{"normal"}: \code{mean}, \code{sd}
#'   \item \code{"uniform"}: \code{min}, \code{max}
#' }
#' They can be given by name or in this order. Draws outside \code{[lower, upper]} are
#' rejected and drawn again (see \code{\link{rsim_prior}}).
#'
#' @param family one of \code{"exponential"}, \code{"lognormal"}, \code{"normal"}
#'   and \code{"uniform"}
#' @param ... the parameters of the family (see Details)
#' @param lower,upper bounds of the draws (default: no truncation)
#'
#' @return an object of class \code{sim_prior}
#' @examples
#' sim_prior("exponential", rate = 2)
#' sim_prior("lognormal", meanlog = -2, sdlog = 1, lower = 1e-8)
#' @export
sim_prior <- function(family = c("exponential", "lognormal", "normal", "uniform"), ...,
                       lower = -Inf, upper = Inf) {
  family <- match.arg(family)
  args <- list(...)
  needed <- switch(family, exponential = "rate", lognormal = c("meanlog", "sdlog"),
                   normal = c("mean", "sd"), uniform = c("min", "max"))
  if (is.null(names(args))) names(args) <- rep("", length(args))
  if (all(names(args) == "") && length(args) == length(needed)) names(args) <- needed
  if (!setequal(names(args), needed) || length(args) != length(needed)) {
    stop("The ", family, " prior needs ", paste0("`", needed, "`", collapse = " and "), ".", call. = FALSE)
  }
  args <- args[needed]
  if (!all(vapply(args, function(a) is.numeric(a) && length(a) == 1 && is.finite(a), logical(1)))) {
    stop("The ", family, " prior needs single finite numbers for ", paste0("`", needed, "`", collapse = " and "), ".", call. = FALSE)
  }
  if (family == "exponential" && args$rate <= 0) stop("The exponential rate must be positive.", call. = FALSE)
  if (family %in% c("lognormal", "normal") && args[[2]] <= 0) stop("The standard deviation of the prior must be positive.", call. = FALSE)
  if (family == "uniform" && args$min >= args$max) stop("The uniform minimum must be smaller than the maximum.", call. = FALSE)
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1 || length(upper) != 1 ||
      is.na(lower) || is.na(upper) || lower >= upper) {
    stop("`lower` must be smaller than `upper`.", call. = FALSE)
  }
  structure(list(family = family, params = args, lower = lower, upper = upper), class = "sim_prior")
}

#' @export
print.sim_prior <- function(x, ...) {
  cat(x$family, "(", paste(names(x$params), "=", format(unlist(x$params)), collapse = ", "), ")", sep = "")
  if (is.finite(x$lower) || is.finite(x$upper)) cat(", truncated to [", x$lower, ", ", x$upper, "]", sep = "")
  cat("\n")
  invisible(x)
}

# the smallest value that a prior can give: the lower end of the support of its family, or its lower bound if that is higher
.support_lower <- function(prior) {
  family_lower <- switch(prior$family, exponential = 0, lognormal = 0, normal = -Inf, uniform = prior$params$min)
  max(prior$lower, family_lower)
}

#' Draw from a prior
#'
#' Draws a vector from a \code{\link{sim_prior}}: the family, its parameters and the
#' bounds say how. The draws outside the bounds are rejected and drawn again
#' (rejection sampling). The function does not set the random seed.
#'
#' @param n number of draws
#' @param prior a \code{\link{sim_prior}}
#' @param max_iter number of rounds of redrawing after which the function gives
#'   up, because the prior and its bounds are (almost) incompatible
#'
#' @return a numeric vector of length \code{n}
#' @examples
#' set.seed(1)
#' rsim_prior(5, sim_prior("normal", mean = 0, sd = 1, lower = 0))
#' @export
rsim_prior <- function(n, prior, max_iter = 1000) {
  stopifnot(inherits(prior, "sim_prior"))
  p <- prior$params
  draw <- function(k) switch(prior$family,
    exponential = stats::rexp(k, p$rate),
    lognormal = stats::rlnorm(k, p$meanlog, p$sdlog),
    normal = stats::rnorm(k, p$mean, p$sd),
    uniform = stats::runif(k, p$min, p$max))
  x <- draw(n)
  bad <- !(x >= prior$lower & x <= prior$upper)
  iter <- 0
  while (any(bad)) {
    iter <- iter + 1
    if (iter > max_iter) {
      stop("The prior and its bounds are (almost) incompatible: no draws fall inside [",
           prior$lower, ", ", prior$upper, "].", call. = FALSE)
    }
    x[bad] <- draw(sum(bad))
    bad <- !(x >= prior$lower & x <= prior$upper)
  }
  x
}

#' Draw probabilities from a Dirichlet distribution
#'
#' Draws a vector of non-negative numbers that sum to 1 from a symmetric Dirichlet
#' distribution, e.g. the frequencies of the outcomes of a barcode edit. With a
#' concentration of 1 every vector of probabilities is equally likely, below 1 the
#' mass is on a few outcomes, above 1 the outcomes are similar. The function does not
#' set the random seed.
#'
#' @param n number of outcomes (components of the vector)
#' @param concentration concentration of the distribution (positive)
#'
#' @return a numeric vector of length \code{n} that sums to 1
#' @examples
#' set.seed(1)
#' sim_dirichlet(3, concentration = 1)
#' @export
sim_dirichlet <- function(n, concentration = 1) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1 || n != round(n)) {
    stop("`n` must be a positive whole number.", call. = FALSE)
  }
  if (!is.numeric(concentration) || length(concentration) != 1 || !(concentration > 0)) {
    stop("`concentration` must be positive.", call. = FALSE)
  }
  for (i in seq_len(100)) {   # a very small concentration can underflow to all zeros
    g <- stats::rgamma(n, concentration)
    if (sum(g) > 0) return(g / sum(g))
  }
  stop("The Dirichlet draw underflowed: use a larger concentration.", call. = FALSE)
}

