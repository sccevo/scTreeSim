#' Prune tree
#' @param obj treedata object 
#' @param rho sampling probability (at present only)
#' @param ntips number of sampled tips
#' @param min_tips minimum number of tips in the pruned tree
#' @param collapse if TRUE, stores the standard structure (node-typed tree), otherwise keeps nodes with one descendant (branch-typed tree)
prune_tree <- function(obj, rho = NA, ntips = NA, min_tips = 2, collapse = TRUE){
  
  # store origin
  origin = obj@phylo$origin
  
  # prune dead particles
  dead_nodes = obj %>% tibble::as_tibble() %>% as.data.frame() %>% dplyr::filter(status == 0) 
  obj = suppressMessages( treeio::drop.tip(obj, dead_nodes$label, collapse.singles = collapse) ) # within the function, there is a tibble command which causes the message "! # Invalid edge matrix for <phylo>. A <tbl_df> is returned." --> suppressed
  
  if (is.null(obj)) {
    message('All particles died! Try another seed.')
    return(NULL)
  }
    
  # prune unsampled particles
  tips = obj@phylo$tip.label
  if (!is.na(rho) & is.na(ntips)) {
    sampled_tips = tips[sample(c(TRUE, FALSE), length(tips), prob = c(rho, 1 - rho), replace = TRUE)]
  } else if (is.na(rho) & !is.na(ntips)) {
    sampled_tips = sample(tips, ntips) 
  } else {
    message('Please provide either the sampling probability or the desired number of tips in pruned tree.')
    return(NULL)
  }

  if (length(sampled_tips) < min_tips) {
    message('Not enough tips sampled! Try another seed.')
    return(NULL)
  }
  
  obj = suppressMessages( treeio::drop.tip(obj, setdiff(tips, sampled_tips), collapse.singles = collapse) )

  # add origin and re-calculate root.edge
  obj@phylo$root.edge = origin - max(ape::node.depth.edgelength(obj@phylo))
  obj@phylo$origin = origin
  
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

