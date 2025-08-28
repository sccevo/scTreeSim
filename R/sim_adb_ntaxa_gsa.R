#' Simulator of a phylogeny under ADB process with GSA
#' 
#' Wrapper for general sampling approach (GSA) for tree simulation with n taxa
#' Compared to the simple sampling approach (SSA), this approach considers subsequent time periods with n particles alive
#' Note: so far, this approach can only be applied to single-type simulations
#' These functions return an object of class phylo
#' 
#' @param ntaxa number of sampled particles (at least 2)
#' @param a scale parameter
#' @param b shape parameter
#' @param d death probability
#' @param rho sampling probability 
#' @param m number of standing taxa that will exist on the first generated trees (before sampling), then be cut; m >> ntaxa/rho
#' @param nsim number of trees to be generated, cut, and sampled from
#' @export
sim_adb_ntaxa_samp_gsa <- function(ntaxa, a, b, d, rho, m = floor(2*ntaxa/rho), nsim = 10) {
  assertthat::assert_that(m > ntaxa)
  
  # generate larger trees
  trees = lapply(c(1:nsim), function(i) {
    suppressMessages({ 
      repeat { 
        tree = sim_adb_ntaxa_complete(ntaxa = m, a = a, b = b, d = d)
        if (!is.null(tree)) break 
      }
      return(tree@phylo) }) })
  
  # apply GSA and prune extinct lineages
  trees_gsa = TreeSim::sim.gsa.taxa(treearray = trees, n = ntaxa/rho, complete = FALSE)
  tree = sample(trees_gsa, 1)[[1]] # sample one tree
  
  # sample tips from tree
  sampled_tips = sample(tree$tip.label, ntaxa) 
  tree = ape::drop.tip(tree, setdiff(tree$tip.label, sampled_tips))
  tree$tip.label = sub('t', '', tree$tip.label)
  
  return(tree)
}


#' Simulator of complete tree under ADB process with GSA
#' @param ntaxa number of sampled particles (at least 2)
#' @param a scale parameter
#' @param b shape parameter
#' @param d death probability
#' @param m number of standing taxa that will exist on the first generated trees, then be cut; m >> ntaxa
#' @param nsim number of trees to be generated, cut, and sampled from
#' @export
sim_adb_ntaxa_complete_gsa <- function(ntaxa, a, b, d, m = 2*ntaxa, nsim = 10) {
  assertthat::assert_that(m > ntaxa)
  
  # generate larger trees
  trees = lapply(c(1:nsim), function(i) {
    suppressMessages({ 
    repeat { 
      tree = sim_adb_ntaxa_complete(ntaxa = m, a = a, b = b, d = d)
      if (!is.null(tree)) break 
    }
    return(tree@phylo) }) })
  
  # apply GSA 
  trees_gsa = TreeSim::sim.gsa.taxa(treearray = trees, n = ntaxa, complete = TRUE)
  tree = sample(trees_gsa, 1)[[1]] # sample one tree
  
  return(tree)
}
  