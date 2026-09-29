#' Simulator of GESTALT-like lineage barcodes
#'
#' Evolves sequences according to the GAPML/ GABI model:
#' Feng, J., DeWitt III, W. S., McKenna, A., Simon, N., Willis, A. D. & Matsen IV, F. A. (2021). Estimation of cell lineage trees by maximum-likelihood phylogenetics. The Annals of Applied Statistics, 15(1), 343 - 362. https://doi.org/10.1214/20-AOAS1400
#' Zwaans, A, Seidel S., Manceau M., Stadler, T. (2025). A Bayesian phylodynamic inference framework for single-cell CRISPR/Cas9 lineage tracing barcode data with dependent target sites. Phil. Trans. R. Soc. B, 380 (1919): 20230318. https://doi.org/10.1098/rstb.2023.0318
#'
#' Cas9 cuts the barcode at one target (focal cut) or at two targets simultaneously (double cut);
#' each cut is repaired by deleting sequence left and/or right of the cut(s) and possibly inserting
#' random nucleotides. An edit that touches a target's cut site deactivates it (it cannot be cut
#' again). Parameters are on the natural scale (Poisson means of the extra length beyond each
#' distribution's minimum).
#'
#' @param tree a treedata object of a sampled tree: ultrametric, i.e. all tips
#'   are sampled cells at present (e.g. from [sim_adb_origin_samp()])
#' @param barcode unedited barcode, a single string of segments (prefix, target_0, spacer,
#'   target_1, ..., suffix) separated by spaces
#' @param clock_rate rate per time unit of edits per target site
#' @param cut_site offset of the cut from the 3' end of each target
#' @param crucial_pos length-2 (left, right) extent of each target's protected region around its cut
#' @param cut_rates per-target cut rate: either one value per target, or a single value used
#'   for every target
#' @param long_trim_factors length-2 (left, right) multiplicative hazard factor for a long trim
#' @param double_cut_weight multiplicative hazard factor for a double (two-target) cut
#' @param insert_zero_prob probability that a cut is repaired with no insertion
#' @param insert_mean Poisson mean of the insertion length beyond 1 (insertion length is
#'   `1 + Poisson(insert_mean)`)
#' @param trim_zero_probs 2x2 matrix (rows left/right, columns focal/double) of the probability
#'   that a cut is repaired with no short deletion on that side (for long deletions this probability is not consulted)
#' @param trim_short_mean length-2 (left, right) Poisson mean of a short trim's length beyond 1
#' @param trim_long_mean length-2 (left, right) Poisson mean of a long trim's length beyond its
#'   long-trim minimum
#'
#' @return a data frame with one row per node in the tree: \code{node}, \code{sequence} (the
#'   edited barcode: original bases uppercase, deleted bases `-`, inserted bases lowercase,
#'   segments space-separated as in `barcode`)
#' @export
sim_barcode_gestalt <- function(tree, barcode, clock_rate = 1, 
                                cut_site = 6, crucial_pos = c(6, 6), cut_rates = 1,
                                long_trim_factors = c(0.1, 0.1), double_cut_weight = 0.1,
                                insert_zero_prob = 0.5, insert_mean = 1,
                                trim_zero_probs = matrix(c(0.1, 0.5, 0.1, 0.3), 2),
                                trim_short_mean = c(1, 3), trim_long_mean = c(8, 4)) {
  
  stopifnot(methods::is(tree, "treedata"))
  if (!ape::is.ultrametric(tree@phylo)) {
    stop("`tree` must be ultrametric (a sampled tree with all tips at present).", call. = FALSE)
  }
  
  barcode <- strsplit(barcode, " ", fixed = TRUE)[[1]]
  layout <- .gestalt_layout(barcode, cut_site, crucial_pos)
  
  if (length(cut_rates) == 1) cut_rates <- rep(cut_rates, layout$n)
  if (length(cut_rates) != layout$n) {
    stop("`cut_rates` must have one entry per target (", layout$n, "), or a single shared value.", call. = FALSE)
  }
  stopifnot(
    "`crucial_pos` must have length 2 (left, right)" = length(crucial_pos) == 2,
    "`long_trim_factors` must have length 2 (left, right)" = length(long_trim_factors) == 2,
    "`trim_zero_probs` must be a 2x2 matrix (rows left/right, columns focal/double)" =
      identical(dim(trim_zero_probs), c(2L, 2L)),
    "`trim_short_mean` must have length 2 (left, right)" = length(trim_short_mean) == 2,
    "`trim_long_mean` must have length 2 (left, right)" = length(trim_long_mean) == 2
  )

  tree_df <- tree %>% tibble::as_tibble() %>% as.data.frame()
  root <- tree_df$node[tree_df$parent == tree_df$node]
  stopifnot(length(root) == 1)
  
  depth_from_root <- ape::node.depth.edgelength(tree@phylo)
  heights <- max(depth_from_root) - depth_from_root
  order_df <- tree_df[order(-heights[tree_df$node]), ]
  
  root_edge <- tree@phylo$root.edge
  origin_height <- if (!is.null(root_edge)) heights[root] + root_edge else heights[root]
  
  evolve <- function(state, branch_length) {
    .gestalt_evolve_branch(state, branch_length, clock_rate, layout, cut_rates, long_trim_factors,
                           double_cut_weight, insert_zero_prob, insert_mean,
                           trim_zero_probs, trim_short_mean, trim_long_mean)
  }
  
  state_at <- list()
  state_at[[as.character(root)]] <- list(
    active = 0:(layout$n - 1),
    allele = list(deleted = rep(FALSE, layout$length), insertions = list())
  )
  
  # evolve along the root/origin edge first, treated like any other branch
  if (!is.null(root_edge) && root_edge > 0) {
    state_at[[as.character(root)]] <- evolve(state_at[[as.character(root)]], origin_height - heights[root])
  }
  
  # then walk every remaining branch, parent state -> child state
  for (i in seq_len(nrow(order_df))) {
    node <- order_df$node[i]
    if (node == root) next
    parent <- order_df$parent[i]
    branch_length <- heights[parent] - heights[node]
    
    state_at[[as.character(node)]] <- evolve(state_at[[as.character(parent)]], branch_length)
  }
  
  sequence <- vapply(
    tree_df$node,
    function(nd) .gestalt_render(state_at[[as.character(nd)]]$allele, barcode),
    character(1)
  )
  
  data.frame(node = tree_df$node, sequence = sequence, stringsAsFactors = FALSE)
}


# Map the cut sites, short and long trim boundaries to positions on the barcode.
# Indexing convention throughout: target indices and barcode positions are 0-based.
# The input barcode here is already separated into segments.
.gestalt_layout <- function(barcode, cut_site, crucial_pos) {
  seg_len <- nchar(barcode)
  n <- (length(barcode) - 1) %/% 2 # number of targets
  total_length <- sum(seg_len)

  # cumulative length up to and including target t
  cum_len <- cumsum(seg_len)
  cut_pos <- cum_len[2 * seq_len(n)] - cut_site

  # protected regions
  prot_left <- pmax(0, cut_pos - crucial_pos[1])
  prot_right <- pmin(total_length - 1, cut_pos + crucial_pos[2] - 1)

  # trim boundaries
  # first/ last targets have no neighbor to reach with a long trim, 
  # so their long-trim minimum is max trim + 1
  Lmax <- c(cut_pos[1], diff(cut_pos) - 1)
  Llong <- c(Lmax[1] + 1, cut_pos[-1] - prot_right[-n])
  Rmax <- c(diff(cut_pos) - 1, total_length - cut_pos[n] - 1)
  Rlong <- c(prot_left[-1] - cut_pos[-n] + 1, Rmax[n] + 1)

  list(
    n = n,
    length = total_length,
    targets = data.frame(
      target = 0:(n - 1),
      cut_pos = cut_pos,
      prot_left = prot_left,
      prot_right = prot_right,
      Lmax = Lmax,
      Llong = Llong,
      Rmax = Rmax,
      Rlong = Rlong
    )
  )
}


# List all possible target tracts given `active` (0-based indices of the currently active targets) out of `n_targets`. 
# A tract cuts at min_cut..max_cut (both active) and deactivates min_deac..max_deac; 
# a long trim on either side is possible regardless of whether the neighbouring target is itself still active 
# (only the barcode ends rule it out).
.gestalt_tracts <- function(active, n_targets) {
  active <- sort(unique(active))

  starts <- do.call(rbind, lapply(active, function(i) {
    opts <- data.frame(min_deac = i, min_cut = i)
    if (i > 0) opts <- rbind(opts, data.frame(min_deac = i - 1, min_cut = i))
    opts
  }))
  ends <- do.call(rbind, lapply(active, function(j) {
    opts <- data.frame(max_cut = j, max_deac = j)
    if (j < n_targets - 1) opts <- rbind(opts, data.frame(max_cut = j, max_deac = j + 1))
    opts
  }))

  tracts <- merge(starts, ends, by = NULL)
  tracts <- tracts[tracts$min_cut <= tracts$max_cut, ]

  tracts$left_long <- tracts$min_deac != tracts$min_cut
  tracts$right_long <- tracts$max_deac != tracts$max_cut
  tracts$focal <- tracts$min_cut == tracts$max_cut

  tracts <- tracts[order(tracts$min_cut, tracts$max_cut, tracts$min_deac, tracts$max_deac), ]
  rownames(tracts) <- NULL
  tracts[, c("min_deac", "min_cut", "max_cut", "max_deac", "left_long", "right_long", "focal")]
}


# Compute the hazard of each tract, given per-target cut rates, (left, right) long-trim factors and the double-cut weight.
# Hazards depend only on the tract (not on the current target status),
# so this could be precomputed once for the fully active status and reused.
.gestalt_hazards <- function(tracts, cut_rates, long_trim_factors, double_cut_weight) {
  side_factor <- long_trim_factors[1]^tracts$left_long * long_trim_factors[2]^tracts$right_long
  lam_min <- cut_rates[tracts$min_cut + 1]
  lam_max <- cut_rates[tracts$max_cut + 1]

  ifelse(tracts$focal,
         lam_min * side_factor,
         (lam_min + lam_max) * double_cut_weight * side_factor)
}


# Draw trim length from bounded shifted Poisson distribution, i.e.
# shift + Poisson(mean), redraw until it does not exceed `max`
.gestalt_shifted_pois <- function(shift, mean, max) {
  repeat {
    len <- shift + stats::rpois(1, mean)
    if (len <= max) return(len)
  }
}

# Draw insertion length from 1 + Poisson(mean) 
# lowercase bases, uniform over a/c/g/t
.gestalt_insertion <- function(mean) {
  len <- 1 + stats::rpois(1, mean)
  paste(sample(c("a", "c", "g", "t"), len, replace = TRUE), collapse = "")
}


# Draw the indel for one chosen `tract` (a single row of .gestalt_tracts()'s output). 
# A long side (left_long/right_long) always gets a deletion, drawn from the long-trim distribution; 
# a short side is optional (trim_zero_probs) and, if it occurs, is drawn from the short-trim distribution, 
# bounded below the neighbouring target's long-trim minimum so it can never reach it. 
# A focal cut with neither side long must leave some trace (redrawn until at least one of insertion/left/right occurs); 
# a double cut always deletes the span between its two cuts regardless, so no such redraw is needed there. 
.gestalt_repair <- function(tract, layout, insert_zero_prob, insert_mean,
                             trim_zero_probs, trim_short_mean, trim_long_mean) {
  min_row <- layout$targets[tract$min_cut + 1, ]
  max_row <- layout$targets[tract$max_cut + 1, ]
  col <- if (tract$focal) 1 else 2

  repeat {
    left_yes <- tract$left_long || (stats::runif(1) > trim_zero_probs[1, col])
    right_yes <- tract$right_long || (stats::runif(1) > trim_zero_probs[2, col])
    insert_yes <- stats::runif(1) > insert_zero_prob
    if (!tract$focal || left_yes || right_yes || insert_yes) break
  }

  left_del <- if (!left_yes) {
    0L
  } else if (tract$left_long) {
    .gestalt_shifted_pois(min_row$Llong, trim_long_mean[1], min_row$Lmax)
  } else {
    .gestalt_shifted_pois(1L, trim_short_mean[1], min_row$Llong - 1)
  }
  right_del <- if (!right_yes) {
    0L
  } else if (tract$right_long) {
    .gestalt_shifted_pois(max_row$Rlong, trim_long_mean[2], max_row$Rmax)
  } else {
    .gestalt_shifted_pois(1L, trim_short_mean[2], max_row$Rlong - 1)
  }
  insert <- if (insert_yes) .gestalt_insertion(insert_mean) else ""

  data.frame(left_del = left_del, right_del = right_del, insert = insert, stringsAsFactors = FALSE)
}


# Apply one drawn `indel` for `tract` to `allele`, i.e. the current state of the sequence 
# (list(deleted = logical over original 0-based positions, insertions = list keyed by the 0-based left-cut position)), 
# returning the updated allele. 
# The deleted span always runs from the left cut minus `left_del` to the right cut plus `right_del` - 1; 
# for a double cut this necessarily includes everything strictly between the two cuts,
# which also removes any earlier insertion anchored in that span 
# (a trim alone never reaches a neighbouring cut, so this only ever triggers via
# a double cut engulfing an already-cut target between its two cuts).
.gestalt_apply <- function(allele, tract, indel, layout) {
  c_min <- layout$targets$cut_pos[tract$min_cut + 1]
  c_max <- layout$targets$cut_pos[tract$max_cut + 1]
  start <- c_min - indel$left_del
  end <- c_max + indel$right_del - 1

  deleted <- allele$deleted
  insertions <- allele$insertions

  if (start <= end) {
    deleted[(start + 1):(end + 1)] <- TRUE
    insertions[names(insertions) %in% as.character(start:end)] <- NULL
  }
  if (nzchar(indel$insert)) insertions[[as.character(c_min)]] <- indel$insert

  list(deleted = deleted, insertions = insertions)
}


# Evolve one lineage's`state` (list(active = 0-based indices of active targets, allele)) 
# along a branch of length `branch_length` by Gillespie simulation. 
# Tracts and their hazards are recomputed from the current `active` set at each iteration.
.gestalt_evolve_branch <- function(state, branch_length, clock_rate, layout, 
                                   cut_rates, long_trim_factors, double_cut_weight, 
                                   insert_zero_prob, insert_mean, 
                                   trim_zero_probs, trim_short_mean, trim_long_mean) {
  active <- state$active
  allele <- state$allele
  t_remaining <- branch_length

  while (length(active) > 0) {
    tracts <- .gestalt_tracts(active, layout$n)
    hazards <- .gestalt_hazards(tracts, cut_rates, long_trim_factors, double_cut_weight)
    total_hazard <- clock_rate * sum(hazards)
    if (total_hazard <= 0) break

    tau <- stats::rexp(1, total_hazard)
    if (tau > t_remaining) break
    t_remaining <- t_remaining - tau

    tract <- tracts[sample.int(nrow(tracts), 1, prob = hazards), ]
    indel <- .gestalt_repair(tract, layout, insert_zero_prob, insert_mean,
                             trim_zero_probs, trim_short_mean, trim_long_mean)
    allele <- .gestalt_apply(allele, tract, indel, layout)
    active <- setdiff(active, tract$min_deac:tract$max_deac)
  }

  list(active = active, allele = allele)
}


# Render an `allele` to the sequence string:
# original bases uppercase, deleted bases "-", inserted bases lowercase, 
# segments (as in `barcode`) separated by spaces. 
# An insertion is placed immediately before the base at its anchor position (the left cut).
.gestalt_render <- function(allele, barcode) {
  bases <- unlist(strsplit(barcode, "", fixed = TRUE))
  seg_end <- cumsum(nchar(barcode))
  seg_start <- seg_end - nchar(barcode) + 1
  chars <- ifelse(allele$deleted, "-", bases)

  segments <- vapply(seq_along(barcode), function(s) {
    idx <- seg_start[s]:seg_end[s]
    pieces <- vapply(idx, function(i) {
      insert <- allele$insertions[[as.character(i - 1)]]
      paste0(if (is.null(insert)) "" else insert, chars[i])
    }, character(1))
    paste(pieces, collapse = "")
  }, character(1))

  paste(segments, collapse = " ")
}

