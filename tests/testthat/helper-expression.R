# a small tree with two types, shared by the expression and counts tests
make_tree <- function() {
  set.seed(42)
  tree <- sim_adb_ntaxa_samp(ntaxa = 8, shape = c(2, 2), mean_lifetime = c(1, 1),
                             death_prob = c(0, 0),
                             asym_trans_prob = matrix(c(0, 0, 0.3, 0), 2, 2),
                             sym_trans_prob = matrix(c(0.7, 0, 0, 1), 2, 2))
  expect_false(is.null(tree))
  tree
}

theta_shared <- function(tree) {
  data.frame(type = sort(unique(tree@data$type)), theta = 0)
}
