
#' Simulator of a phylogeny from an Age-Dependent Branching Process for a fixed number of sampled particles
#' @param ntaxa number of sampled particles (at least 2)
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param rho sampling probability 
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param Xi_as matrix of asymetric type transition probabilities
#' @param Xi_s matrix of symetric type transition probabilities
#' @export
sim_adb_ntaxa_samp <- function(ntaxa, a, b, d = 0, rho = 1, origin_type = 0, Xi_as = matrix(0), Xi_s = matrix(1)) {
  # assert that all inputs are correct
  ntypes = length(a)
  assert_that(all(c(length(b) == ntypes, length(d) == ntypes, 
                    dim(Xi_as) == c(ntypes, ntypes), dim(Xi_s) == c(ntypes, ntypes),
                    ntaxa >= 2, origin_type %in% c(0:ntypes),
                    d >= 0, d < 1, rho > 0, rho <= 1)),
              msg = 'The inputs do not have proper dimensions or values. Please check all parameters!')
  assert_that(all(sapply(seq(1, ntypes), function(i) {sum(2*Xi_as[i, ]) + sum(Xi_s[i, ]) == 1})),
              msg = 'The transition probabilities do not some to 1. Please check!')
  
  # estimate the number of taxa in the full tree
  nfull = ntaxa / rho
  
  # simulate full tree
  tree = sim_adb_ntaxa_complete(nfull, a, b, d, origin_type, Xi_as, Xi_s)
  if (is.null(tree)) {
    return(NULL)
  }
  
  phylogeny <- prune_tree(tree, rho, min_tips)

  return(phylogeny)
}


#' Simulator of the complete Age-Dependent Branching Process (up to a fixed number of living particles)
#' @param ntaxa number of sampled particles (at least 2)
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param Xi_as matrix of asymetric type transition probabilities
#' @param Xi_s matrix of symetric type transition probabilities
#' @export
sim_adb_ntaxa_complete <- function(ntaxa, a, b, d, origin_type=0, Xi_as=matrix(0), Xi_s=matrix(0)) {
  
  # initialize
  edges = matrix(nrow = 0, ncol = 2)
  edge_lengths = c()
  nodes = data.frame(
    id = integer(0), 
    height = numeric(0), 
    type = integer(0), 
    parent = integer(0), 
    leftchild = integer(0), 
    rightchild = integer(0), 
    status = integer(0) # 0 = dead, 1 = alive, 2 = divided
  )
  
  # sample the lifetime of the first particle
  root_edge = rgamma(1, shape = b[origin_type + 1], scale = a[origin_type + 1])
  nodes = bind_rows(nodes, c(id = 1, height = root_edge, type = origin_type, parent = NA, leftchild = NA, rightchild = NA, status = 1))
  events = nodes
  event_counter = 1
  living = 1
  
  while (living < ntaxa & nrow(events) > 0) {
    if (living %% 10000 == 0) {
      print(paste0(living, ' cells alive'))
    }
    event = as.list(events[1, ]); events = events[-1, ] # look at one event
    
    if (runif(1) < d[event$type + 1]) {
      # particle dies
      nodes[event$id, "status"] = 0
      living = living - 1
    } else {
      # particle divides, create two new children
      left_id = event_counter + 1
      right_id = event_counter + 2
      event_counter = event_counter + 2
      nodes[event$id, "status"] = 2
      living = living + 1

      if (ncol(Xi_s) == 1) { 
        # single-type case
        children_types = rep(origin_type, 2)
      } else {
        # multi-type case: sample types
        children_types = sample_types(event$type, Xi_as, Xi_s)
      }
      
      # sample lifetimes and add new nodes
      left_type = children_types[1]
      left_lifetime = rgamma(1, shape = b[left_type + 1], scale = a[left_type + 1])
      left_node = c(id = left_id, height = event$height + left_lifetime, type = left_type, 
                    parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      nodes = bind_rows(nodes, left_node)

      right_type = children_types[2]
      right_lifetime = rgamma(1, shape = b[right_type + 1], scale = a[right_type + 1])
      right_node = c(id = right_id, height = event$height + right_lifetime, type = right_type, 
                     parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      nodes = bind_rows(nodes, right_node)
      
      # to be processed, sort by height
      events = bind_rows(events, left_node, right_node) 
      events = events %>% arrange(height)
      
      # add child relationships and append the current edges and edge lengths
      nodes[event$id, c("leftchild", "rightchild")] = c(left_id, right_id)
      edges = rbind(edges, c(event$id, left_id), c(event$id, right_id))
    }
  }
  
  if (sum(nodes$status == 1) < ntaxa) {
    message('Too many particles died! Try another seed or decrease death probability.')
    return(NULL)
  }
  
  # truncate to the minimal height of living particles
  stopping_time = min(nodes[nodes$status == 1, 'height'])
  assert_that(stopping_time > max(nodes[nodes$status == 2, 'height'])) # ensure that truncation occurs after any cell division
  nodes[nodes$status == 1, 'height'] = stopping_time # update heights
  
  # calculate edge lengths post-hoc
  edge_lengths = sapply(c(1:nrow(edges)), function(i) {
    nodes[nodes$id == edges[i, 2], 'height'] - nodes[nodes$id == edges[i, 1], 'height'] })
  
  # check number of tips
  Nnode = sum(nodes$status == 2)
  Ntip = sum(nodes$status %in% c(0,1))
  
  # assign labels
  nodes[which(nodes$status %in% c(0,1)), 'label'] = c(1:Ntip) 
  nodes[which(nodes$status == 2), 'label'] = c((Ntip + 1):(Ntip + Nnode)) 
  
  # use labels in edge matrix
  edges_recoded = matrix(nodes[as.double(edges), 'label'], nrow = nrow(edges), ncol = ncol(edges))
  
  # create phylogenetic tree
  phylo_tree = list(edge = edges_recoded, edge.length = edge_lengths, Nnode = Nnode, 
                    tip.label = as.character(nodes[which(nodes$status %in% c(0,1)), 'label']))
  class(phylo_tree) = "phylo"
  
  # create treedata object
  tree = as.treedata(phylo_tree)
  tree@phylo$root.edge = root_edge 
  tree@phylo$origin = stopping_time
  data = as_tibble(tree)
  types = nodes %>% 
    select(node = label, status, type) %>% 
    mutate(type = as.factor(type)) %>%
    arrange(node) %>%
    as_tibble
  tree@data = types
  
  return(tree)
}

  
