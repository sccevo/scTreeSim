#' Prune tree
#' @param obj S4 object with @phylo class 
#' @param rho sampling probability (at present only)
prune_tree <- function(obj, rho, min_tips){

  # prune dead particles
  dead_nodes = obj %>% as_tibble %>% filter(status == 0) %>% pull(label)
  phylogeny = drop.tip(obj@phylo, dead_nodes)
  if (is.null(phylogeny)) {
    message('All particles died! Try another seed.')
    return(NULL)
  }

  # prune unsampled particles
  tips = phylogeny$tip.label
  sampled_tips = sample(c(TRUE, FALSE), length(tips), prob = c(rho, 1 - rho), replace = TRUE)
  if (sum(sampled_tips) < min_tips) {
    message('Not enough tips sampled! Try another seed.')
    return(NULL)
  }
  phylogeny = drop.tip(phylogeny, tips[!sampled_tips])
#  phylogeny$tip.label = as.character(c(1:length(phylogeny@phylo$tip.label)))
  
  # remove height from data for exporting tree
  #phylogeny@data = phylogeny@data %>% select(-status)
                                             
  # keep origin in final object
 # phylogeny@phylo$origin = tree@phylo$origin
  obj@phylo = phylogeny
  return(obj)

}



#' Helper function for sampling types of offspring upon division
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

