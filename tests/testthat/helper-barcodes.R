# Helpers for comparing the barcode simulators on a fixed small tree.
#
# Tree ((a:1, b:1):0.5, c:1.5) with a root edge of 0.5: every tip is at distance
# 2 from the origin. Nodes: tips 1-3, root 4, internal node 5.

small_barcode_tree <- function() {
  phy <- ape::read.tree(text = "((a:1,b:1):0.5,c:1.5);")
  phy$root.edge <- 0.5
  tr <- treeio::as.treedata(phy)
  tr@data <- tibble::tibble(node = 1:5, status = c(1, 1, 1, 2, 2), type = 0L)
  tr
}
small_barcode_time <- 2  # time from the origin to the tips

# Closed-form distribution of a single site at a tip after time `t`:
# unedited, each edit outcome, silenced (edit and silencing are independent clocks).
single_site_probs <- function(edit_rate, edit_probs, silencing_rate, t = small_barcode_time) {
  c("0" = exp(-(edit_rate + silencing_rate) * t),
    setNames(edit_probs * (exp(-silencing_rate * t) - exp(-(edit_rate + silencing_rate) * t)),
             seq_along(edit_probs)),
    "-" = 1 - exp(-silencing_rate * t))
}

# Single-site states (strings) of every node, for each replicate (column), from any of the
# three simulators; the replicates are independent barcodes / targets of one simulator call.
replicate_states <- function(out) {
  m <- as.matrix(out[-1])
  storage.mode(m) <- "character"
  m
}

# Count tuples of states of the given rows (nodes) over the replicates
joint_counts <- function(states, nodes) {
  table(do.call(paste, c(lapply(nodes, function(nd) states[nd, ]), sep = "|")))
}

# Chi-square test of homogeneity between two simulators' joint counts
# (Monte Carlo p-value: small cells are common)
homogeneity_p <- function(counts_a, counts_b) {
  cells <- union(names(counts_a), names(counts_b))
  tab <- rbind(counts_a[cells], counts_b[cells])
  tab[is.na(tab)] <- 0
  suppressWarnings(stats::chisq.test(tab, simulate.p.value = TRUE, B = 5000)$p.value)
}

# Four binomial standard errors
tolerance_4se <- function(p, n) 4 * sqrt(p * (1 - p) / n)

# Observed frequencies are within `tol` (absolute) of the expected ones
expect_within <- function(observed, expected, tol) {
  testthat::expect_true(all(abs(as.numeric(observed) - as.numeric(expected)) <= tol),
                        info = paste("max deviation", max(abs(as.numeric(observed) - as.numeric(expected))),
                                     "> tolerance", max(tol)))
}
