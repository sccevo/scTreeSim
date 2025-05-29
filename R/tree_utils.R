#' Prune tree
#' @param obj treedata object 
#' @param rho sampling probability (at present only)
#' @param min_tips minimum number of tips in the pruned tree
prune_tree <- function(obj, rho, min_tips){
  # get root edge
  stem = obj@phylo$root.edge
  # prune dead particles
  dead_nodes = obj %>% tibble::as_tibble() %>% dplyr::filter(status == 0) %>% dplyr::pull(label)
  obj = treeio::drop.tip(obj, dead_nodes)

  if (is.null(obj@phylo)) {
    message('All particles died! Try another seed.')
    return(NULL)
  }
    
  # prune unsampled particles
  tips = obj@phylo$tip.label
  sampled_tips = sample(c(TRUE, FALSE), length(tips), prob = c(rho, 1 - rho), replace = TRUE)
  if (sum(sampled_tips) < min_tips) {
    message('Not enough tips sampled! Try another seed.')
    return(NULL)
  }
  obj = treeio::drop.tip(obj, tips[!sampled_tips])
 
  # add back root edge
  obj@phylo$root.edge = stem
  return(obj)

}



#' Helper function for sampling types of offspring upon division
#' @param parent_type ID of parent type
#' @param Xi_as matrix of asymmetric type transition probabilities
#' @param Xi_s matrix of symmetric type transition probabilities
sample_types <- function(parent_type, Xi_as, Xi_s) {
  ntypes = ncol(Xi_s)
  
  cum_prob = 0
  for (i in seq(0, ntypes - 1)) {
    # symmetric case
    cum_prob = cum_prob + Xi_s[parent_type + 1, i + 1];
    if (runif(1) < cum_prob) {
      return(c(i, i))
    }
    # asymmetric case
    cum_prob = cum_prob + 2*Xi_as[parent_type + 1, i + 1];
    if (runif(1) < cum_prob) {
      if (runif(1) < 0.5) {
        return(c(parent_type, i))
      } else {
        return(c(i, parent_type))
      }
    }
  }
}

