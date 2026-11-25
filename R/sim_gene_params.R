# the root type and the height of a tree
.tree_facts <- function(tree) {
  phy <- tree@phylo
  root <- setdiff(phy$edge[, 1], phy$edge[, 2])
  list(root_type = tree@data$type[tree@data$node == root], height = max(ape::node.depth.edgelength(phy)))
}

# the number of types of a tree, unless it is given
.n_types <- function(tree, n_types) {
  if (is.null(n_types)) n_types <- max(tree@data$type) + 1L
  if (!is.numeric(n_types) || length(n_types) != 1 || is.na(n_types) || n_types < 1 || n_types != round(n_types)) {
    stop("`n_types` must be a positive whole number.", call. = FALSE)
  }
  if (!all(tree@data$type %in% (seq_len(n_types) - 1L))) {
    stop("The tree has types outside 0, ..., ", n_types - 1, ".", call. = FALSE)
  }
  as.integer(n_types)
}

# The data frame of true parameters of the genes, shared by sim_gene_params() and gene_params_from_table()
.gene_truth_table <- function(gene, gene_class, alpha, sigma2, sd_stat, theta, x_root, theta_base, log_g,
                              eta = NULL, psi = NULL, gap_cutoff = 1) {
  n_types <- ncol(theta)
  theta_gap <- apply(theta, 1, max) - apply(theta, 1, min)
  class_eff <- ifelse(gene_class == "BM", "BM", ifelse(theta_gap > gap_cutoff * sd_stat, "OUM", "OU1"))
  truth <- data.frame(gene = gene, class = gene_class, class_eff = class_eff, theta_gap = theta_gap,
                      log_g = log_g, alpha = alpha, sigma2 = sigma2, sd_stat = sd_stat, theta_base = theta_base)
  if (!is.null(eta)) truth$eta <- eta
  if (!is.null(psi)) for (k in seq_len(n_types - 1)) truth[[paste0("psi_", k - 1)]] <- psi[, k]
  for (k in seq_len(n_types)) truth[[paste0("theta_", k - 1)]] <- theta[, k]
  truth$x_root <- x_root
  truth
}

#' Draw the parameters of many genes from priors
#'
#' Draws the parameters of an independent BM / OU1 / OUM expression process for
#' every gene, to be simulated on a tree with \code{\link{sim_expression_matrix}}.
#' The function does not set the random seed: call \code{set.seed()} before it
#' to get the same genes again.
#'
#' Every gene is an Ornstein-Uhlenbeck (OU) process with probability \code{pi_ou}
#' and a Brownian motion (BM, \code{alpha = 0}) otherwise. The regimes are the
#' types \code{0, ..., n_types - 1} of the tree, and the optimum of type \code{r} is
#' \deqn{\theta_r = \theta_{base} + \eta (R \psi_r - 1),}
#' where \eqn{R} is the number of types, \eqn{\psi} is drawn from a symmetric
#' Dirichlet distribution with concentration \code{psi_concentration} and
#' \eqn{\theta_{base} = \log g}, with \eqn{\log g} normal with the mean that gives
#' \code{mean_count_gene} as the mean of \eqn{g} and standard deviation \code{log_g_sd}
#' (the baseline level of the gene). \eqn{\eta = 0} gives the same optimum to
#' every type (OU1); larger \eqn{\eta} separates them (OUM). The root value of an
#' OU gene is drawn from the stationary distribution of the root's type, and that
#' of a BM gene is \eqn{\theta_{base}}.
#'
#' @param tree a treedata object, e.g. from \code{\link{sim_adb_ntaxa_samp}}, with
#'   the types \code{0, ..., n_types - 1}
#' @param n_genes number of genes
#' @param n_types number of types (regimes). Default: the largest type of the tree
#'   plus one.
#' @param pi_ou probability that a gene is OU rather than BM
#' @param alpha,sigma2,eta \code{\link{sim_prior}}s of the strength of selection
#'   (OU genes only), the variance of the driving noise and the separation of the
#'   optima (OU genes only). None of them can be negative: a prior that can give
#'   negative values (a normal, say) needs a lower bound of at least 0. The default bounds are those of the app: an OU gene
#'   has \code{alpha >= 0.1} and the variance is at most 2.
#' @param psi_concentration concentration of the Dirichlet distribution of the
#'   share of each type in the deviations of the optima: below 1 the optima are far
#'   apart, above 1 close
#' @param mean_count_gene mean of the baseline level \eqn{g} of the genes
#' @param log_g_sd standard deviation of \eqn{\log g} across genes
#' @param gap_cutoff an OU gene is effectively OUM if the largest difference between
#'   its optima is more than \code{gap_cutoff} times its stationary standard
#'   deviation, and OU1 otherwise (descriptive only)
#'
#' @return A data frame with one row per gene and the columns \code{gene} (names
#'   \code{gene_1}, ...), \code{class} (BM or OU), \code{class_eff} (BM, OU1 or OUM),
#'   \code{theta_gap}, \code{log_g}, \code{alpha}, \code{sigma2}, \code{sd_stat}
#'   (the stationary standard deviation of an OU gene, and the standard deviation at
#'   the tips of a BM gene), \code{theta_base}, \code{eta}, \code{psi_0}, ... (the
#'   shares of all types but the last), \code{theta_0}, \code{theta_1}, ... (the
#'   optima) and \code{x_root}. Pass \code{alpha}, \code{sigma2}, \code{x_root} and
#'   \code{\link{gene_optima}} of it to \code{\link{sim_expression_matrix}}.
#' @examples
#' set.seed(1)
#' tree <- sim_adb_ntaxa_samp(20, scale = c(0.2, 0.2), shape = c(5, 5), death_prob = c(0.05, 0.05),
#'                            asym_trans_prob = matrix(c(0, 0, 0.2, 0), 2),
#'                            sym_trans_prob = matrix(c(0.8, 0, 0, 1), 2), collapse = FALSE)
#' genes <- sim_gene_params(tree, n_genes = 10)
#' head(genes[c("gene", "class_eff", "alpha", "theta_0", "theta_1")], 3)
#' expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2,
#'                               root_value = genes$x_root, theta = gene_optima(genes),
#'                               gene_names = genes$gene)
#' dim(expr)
#' @export
sim_gene_params <- function(tree, n_genes = 100, n_types = NULL, pi_ou = 0.5,
                            alpha = sim_prior("exponential", rate = 1, lower = 0.1),
                            sigma2 = sim_prior("exponential", rate = 1, lower = 1e-8, upper = 2),
                            eta = sim_prior("exponential", rate = 1, lower = 0),
                            psi_concentration = 0.5, mean_count_gene = 1, log_g_sd = 1.5, gap_cutoff = 1) {
  stopifnot(methods::is(tree, "treedata"))
  if (!is.numeric(n_genes) || length(n_genes) != 1 || is.na(n_genes) || n_genes < 1 || n_genes != round(n_genes)) {
    stop("`n_genes` must be a positive whole number.", call. = FALSE)
  }
  if (!is.numeric(pi_ou) || length(pi_ou) != 1 || is.na(pi_ou) || pi_ou < 0 || pi_ou > 1) {
    stop("`pi_ou` must be a probability.", call. = FALSE)
  }
  for (p in list(alpha, sigma2, eta)) {
    if (!inherits(p, "sim_prior")) stop("`alpha`, `sigma2` and `eta` must be created with sim_prior().", call. = FALSE)
    if (.support_lower(p) < 0) stop("The priors of `alpha`, `sigma2` and `eta` need a lower bound of at least 0: they cannot be negative.", call. = FALSE)
  }
  if (!is.numeric(psi_concentration) || length(psi_concentration) != 1 || !(psi_concentration > 0)) {
    stop("`psi_concentration` must be positive.", call. = FALSE)
  }
  if (!is.numeric(mean_count_gene) || length(mean_count_gene) != 1 || !(mean_count_gene > 0)) {
    stop("`mean_count_gene` must be positive.", call. = FALSE)
  }
  if (!is.numeric(log_g_sd) || length(log_g_sd) != 1 || !(log_g_sd >= 0)) {
    stop("`log_g_sd` must be non-negative.", call. = FALSE)
  }
  n_types <- .n_types(tree, n_types)
  facts <- .tree_facts(tree)
  root_regime <- facts$root_type + 1L
  log_g_mean <- log(mean_count_gene) - log_g_sd^2 / 2   # E[g] = exp(mu + sd^2 / 2) = mean_count_gene

  is_bm <- stats::rbinom(n_genes, 1, pi_ou) == 0
  gene_class <- ifelse(is_bm, "BM", "OU")

  log_g   <- stats::rnorm(n_genes, log_g_mean, log_g_sd)
  alpha   <- ifelse(is_bm, 0, rsim_prior(n_genes, alpha))
  sigma2  <- rsim_prior(n_genes, sigma2)
  # sd of x at the tips (BM) / stationary sd (OU)
  sd_stat <- ifelse(is_bm, sqrt(sigma2 * facts$height), sqrt(sigma2 / (2 * alpha)))

  eta <- ifelse(is_bm, 0, rsim_prior(n_genes, eta))   # BM has no optima: eta and psi are unused
  # psi: the share of each type in the deviations of the optima, Dirichlet(psi_concentration, ...)
  psi <- if (n_types == 1) {
    matrix(1, n_genes, 1)
  } else if (n_types == 2) {
    psi1 <- stats::rbeta(n_genes, psi_concentration, psi_concentration)
    cbind(psi1, 1 - psi1)
  } else {
    g <- matrix(stats::rgamma(n_genes * n_types, psi_concentration), n_genes, n_types)
    g / rowSums(g)
  }
  psi[is_bm, ] <- 1 / n_types
  dev <- eta * (n_types * psi - 1)                     # deviations of the optima (they sum to 0 per gene)

  # theta_base = log g: g is the baseline level of the gene (the median of its noise-free counts)
  theta_base <- log_g
  theta <- theta_base + dev                            # genes x types
  x_root <- ifelse(is_bm, theta_base,
                   stats::rnorm(n_genes, theta[cbind(seq_len(n_genes), root_regime)], sd_stat))

  .gene_truth_table(paste0("gene_", seq_len(n_genes)), gene_class, alpha, sigma2, sd_stat, theta, x_root,
                    theta_base = theta_base, log_g = log_g, eta = eta, psi = psi, gap_cutoff = gap_cutoff)
}

#' Parameters of genes from a table
#'
#' The same table of true parameters as \code{\link{sim_gene_params}}, for genes
#' whose parameters you specify yourself (for example read with
#' \code{read.csv}).
#'
#' @param table a data frame with one row per gene and the columns
#'   \code{alpha} (strength of selection, \code{0} is a Brownian motion),
#'   \code{sigma2} (variance of the driving noise, positive) and
#'   \code{theta_0}, \code{theta_1}, ... (the optimum of each type; they may be
#'   missing for Brownian motions). Optional columns: \code{gene} (names; default
#'   \code{gene_1}, ...) and \code{root_value} (default: the optimum of the type of
#'   the root for OU genes and 0 for Brownian motions).
#' @inheritParams sim_gene_params
#'
#' @return A data frame with the columns that \code{\link{sim_gene_params}}
#'   returns, except \code{eta} and \code{psi_*}. \code{theta_base}
#'   (\code{log_g}) is the mean of the optima of an OU gene and the root value of
#'   a BM gene.
#' @examples
#' set.seed(1)
#' tree <- sim_adb_ntaxa_samp(20, scale = c(0.2, 0.2), shape = c(5, 5), death_prob = c(0.05, 0.05),
#'                            asym_trans_prob = matrix(c(0, 0, 0.2, 0), 2),
#'                            sym_trans_prob = matrix(c(0.8, 0, 0, 1), 2), collapse = FALSE)
#' table <- data.frame(alpha = c(0, 1, 2), sigma2 = 0.5, theta_0 = c(NA, 0, -1), theta_1 = c(NA, 0, 1))
#' gene_params_from_table(table, tree)[c("gene", "class_eff", "theta_gap")]
#' @export
gene_params_from_table <- function(table, tree, n_types = NULL, gap_cutoff = 1) {
  stopifnot(methods::is(tree, "treedata"))
  tab <- as.data.frame(table, stringsAsFactors = FALSE)
  if (nrow(tab) == 0) stop("The table has no rows.", call. = FALSE)
  names(tab) <- trimws(names(tab))
  names(tab)[names(tab) == "sigmasq"] <- "sigma2"
  n_types <- .n_types(tree, n_types)
  theta_cols <- paste0("theta_", seq_len(n_types) - 1L)
  missing_cols <- setdiff(c("alpha", "sigma2"), names(tab))
  if (length(missing_cols) > 0) {
    stop("The table needs the columns alpha and sigma2 (missing: ", paste(missing_cols, collapse = ", "), ").", call. = FALSE)
  }
  for (cn in intersect(c("alpha", "sigma2", "root_value", theta_cols), names(tab))) {
    tab[[cn]] <- suppressWarnings(as.numeric(tab[[cn]]))
  }
  if (anyNA(tab$alpha) || any(tab$alpha < 0)) stop("alpha must be a number >= 0 for every gene (0 is a Brownian motion).", call. = FALSE)
  if (anyNA(tab$sigma2) || any(tab$sigma2 <= 0)) stop("sigma2 must be a number > 0 for every gene.", call. = FALSE)
  is_bm <- tab$alpha == 0
  absent <- setdiff(theta_cols, names(tab))
  if (length(absent) > 0) {
    if (!all(is_bm)) {
      stop("The tree has ", n_types, " type(s): the table needs the columns ", paste(theta_cols, collapse = ", "),
           " (missing: ", paste(absent, collapse = ", "), ").", call. = FALSE)
    }
    for (cn in absent) tab[[cn]] <- 0   # only Brownian motions: the optima are not used
  }
  for (cn in theta_cols) {
    if (any(is.na(tab[[cn]]) & !is_bm)) stop("The optimum ", cn, " is missing for some OU genes (alpha > 0).", call. = FALSE)
    tab[[cn]][is.na(tab[[cn]])] <- 0
  }
  if (!"gene" %in% names(tab)) tab$gene <- paste0("gene_", seq_len(nrow(tab)))
  tab$gene <- as.character(tab$gene)
  if (anyDuplicated(tab$gene) || any(tab$gene == "")) stop("The gene names must be unique and not empty.", call. = FALSE)
  if ("root_value" %in% names(tab) && anyNA(tab$root_value)) {
    stop("root_value has missing values: leave the column out to use the default.", call. = FALSE)
  }

  facts <- .tree_facts(tree)
  theta <- as.matrix(tab[theta_cols])
  dimnames(theta) <- NULL
  x_root <- if ("root_value" %in% names(tab)) tab$root_value else ifelse(is_bm, 0, theta[, facts$root_type + 1L])
  sd_stat <- ifelse(is_bm, sqrt(tab$sigma2 * facts$height), sqrt(tab$sigma2 / (2 * tab$alpha)))
  theta_base <- ifelse(is_bm, x_root, rowMeans(theta))
  .gene_truth_table(tab$gene, ifelse(is_bm, "BM", "OU"), tab$alpha, tab$sigma2, sd_stat, theta, x_root,
                    theta_base = theta_base, log_g = theta_base, gap_cutoff = gap_cutoff)
}

#' Optima of the genes, as the \code{theta} of the expression simulation
#'
#' Arranges the columns \code{theta_0}, \code{theta_1}, ... of the parameters of
#' many genes (\code{\link{sim_gene_params}}, \code{\link{gene_params_from_table}})
#' as the \code{theta} argument of \code{\link{sim_expression_matrix}}.
#'
#' @param gene_params a data frame of parameters of genes, with the columns
#'   \code{gene} and \code{theta_0}, \code{theta_1}, ...
#'
#' @return a data frame with the column \code{type} (\code{0, 1, ...}) and one
#'   column per gene, named by \code{gene_params$gene}, with the optimum of the
#'   gene in each type
#' @export
gene_optima <- function(gene_params) {
  cols <- grep("^theta_[0-9]+$", names(gene_params), value = TRUE)
  if (length(cols) == 0) stop("`gene_params` needs the columns theta_0, theta_1, ...", call. = FALSE)
  cols <- cols[order(as.integer(sub("theta_", "", cols)))]
  theta <- data.frame(type = seq_along(cols) - 1L, t(as.matrix(gene_params[cols])))
  names(theta) <- c("type", gene_params$gene)
  theta
}

#' Draw the parameters of the counts of many genes and cells
#'
#' Draws what \code{\link{sim_expression_counts_matrix}} needs: the dispersion
#' \code{phi} of every gene and the size factor of every cell. The function does
#' not set the random seed: call \code{set.seed()} before it to get the same
#' values again.
#'
#' The size factors are lognormal with mean 1. The dispersion is either drawn from
#' a prior, independently of the expression (\code{phi_mode = "prior"}), or
#' follows a trend that decreases with the baseline level of the gene
#' (\code{phi_mode = "trend"}):
#' \deqn{\phi_g = (a_0 + a_1 / g_g) e^{z_g}, \quad z_g \sim N(0, sd^2).}
#'
#' With \code{target_mean_count}, a constant is added to the expression of every
#' gene so that the expected mean count over all genes and cells is exactly the target.
#' The constant is also added to the log-scale columns of the returned parameters
#' (\code{log_g}, \code{theta_base}, \code{theta_0}, ..., \code{x_root}).
#'
#' @param expr data frame of the expression of the genes on the tree, from
#'   \code{\link{sim_expression_matrix}}
#' @param tree the treedata object of \code{expr}
#' @param gene_params the parameters of the genes of \code{expr}
#'   (\code{\link{sim_gene_params}} or \code{\link{gene_params_from_table}}); the
#'   column \code{gene} gives the genes and the column \code{log_g} is used by the
#'   trend
#' @param cell_size_sd standard deviation (log scale) of the size factors
#' @param target_mean_count mean count per gene and cell to calibrate to, or
#'   \code{NULL} (default) for no calibration
#' @param phi a \code{\link{sim_prior}} of the dispersion (\code{phi_mode = "prior"})
#' @param phi_mode \code{"prior"} or \code{"trend"}
#' @param phi_a0,phi_a1,phi_sd parameters of the trend of the dispersion
#'
#' @return A list with \code{expr} (the expression, shifted when
#'   \code{target_mean_count} is given), \code{phi} (one value per gene),
#'   \code{size_factor} (one value per sampled tip, named by the node id of the
#'   tip), \code{shift} (the constant added to the expression) and
#'   \code{gene_params} (\code{gene_params}, shifted, with the column \code{phi}).
#'   Pass the first three to \code{\link{sim_expression_counts_matrix}}.
#' @examples
#' set.seed(1)
#' tree <- sim_adb_ntaxa_samp(20, scale = c(0.2, 0.2), shape = c(5, 5), death_prob = c(0.05, 0.05),
#'                            asym_trans_prob = matrix(c(0, 0, 0.2, 0), 2),
#'                            sym_trans_prob = matrix(c(0.8, 0, 0, 1), 2), collapse = FALSE)
#' genes <- sim_gene_params(tree, n_genes = 10)
#' expr <- sim_expression_matrix(tree, alpha = genes$alpha, sigmasq = genes$sigma2,
#'                               root_value = genes$x_root, theta = gene_optima(genes),
#'                               gene_names = genes$gene)
#' counts_params <- sim_count_params(expr, tree, genes, target_mean_count = 2)
#' counts <- sim_expression_counts_matrix(counts_params$expr, tree, phi = counts_params$phi,
#'                                        size_factor = counts_params$size_factor)
#' dim(counts)
#' @export
sim_count_params <- function(expr, tree, gene_params, cell_size_sd = 0.9, target_mean_count = NULL,
                             phi = sim_prior("lognormal", meanlog = -2, sdlog = 2, lower = 1e-8),
                             phi_mode = c("prior", "trend"), phi_a0 = 0.1, phi_a1 = 0.5, phi_sd = 0.3) {
  stopifnot(methods::is(tree, "treedata"), is.data.frame(expr), is.data.frame(gene_params))
  phi_mode <- match.arg(phi_mode)
  if (phi_mode == "prior") {
    if (!inherits(phi, "sim_prior")) stop("`phi` must be created with sim_prior().", call. = FALSE)
    if (.support_lower(phi) < 0) stop("The prior of `phi` needs a lower bound of at least 0: the dispersion cannot be negative.", call. = FALSE)
  }
  if (!is.numeric(cell_size_sd) || length(cell_size_sd) != 1 || !(cell_size_sd >= 0)) {
    stop("`cell_size_sd` must be non-negative.", call. = FALSE)
  }
  if (!is.null(target_mean_count) && (!is.numeric(target_mean_count) || length(target_mean_count) != 1 || !(target_mean_count > 0))) {
    stop("`target_mean_count` must be positive.", call. = FALSE)
  }
  gene_ids <- gene_params$gene
  if (!all(gene_ids %in% names(expr))) stop("`expr` needs a column for every gene of `gene_params`.", call. = FALSE)
  tips <- tree@data$node[tree@data$status == 1]            # the sampled tips are the cells
  n_tip <- length(tips)
  n_gene <- length(gene_ids)

  # in "trend" mode phi is computed after the level calibration, as it depends on the final log_g
  phi_draw <- if (phi_mode == "prior") rsim_prior(n_gene, phi) else stats::rnorm(n_gene, 0, phi_sd)
  size_factor <- exp(stats::rnorm(n_tip, 0, cell_size_sd))
  size_factor <- size_factor / mean(size_factor)           # mean 1: g stays on the counts scale
  names(size_factor) <- as.character(tips)

  # exact level calibration: the constant c such that mean(exp(x_tip + c) * size_factor) = target
  shift <- 0
  if (!is.null(target_mean_count)) {
    x_tip <- t(as.matrix(expr[match(tips, expr$node), gene_ids]))   # genes x cells
    shift <- log(target_mean_count / mean(exp(x_tip) * matrix(size_factor, n_gene, n_tip, byrow = TRUE)))
    expr[gene_ids] <- expr[gene_ids] + shift
    shifted <- intersect(c("log_g", "theta_base", "x_root", grep("^theta_[0-9]+$", names(gene_params), value = TRUE)),
                         names(gene_params))
    gene_params[shifted] <- gene_params[shifted] + shift
  }
  gene_params$phi <- if (phi_mode == "prior") {
    phi_draw
  } else {
    if (!"log_g" %in% names(gene_params)) stop("`gene_params` needs the column log_g for the trend.", call. = FALSE)
    (phi_a0 + phi_a1 / exp(gene_params$log_g)) * exp(phi_draw)
  }
  list(expr = expr, phi = gene_params$phi, size_factor = size_factor, shift = shift, gene_params = gene_params)
}
