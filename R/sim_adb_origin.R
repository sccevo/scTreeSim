
#' Simulator of a phylogeny from an Age-Dependent Branching Process for a fixed time interval (since origin)
#' @param origin_time time of birth of the initial particle
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param rho sampling probability 
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param Xi_as matrix of asymetric type transition probabilities
#' @param Xi_s matrix of symetric type transition probabilities
#' @param min_tips minimum number of tips in the phylogeny
#' @export
sim_adb_origin_samp <- function(origin_time, a, b, d = 0, rho = 1, origin_type = 0, Xi_as = matrix(0), Xi_s = matrix(1), min_tips = 2) {
  # assert that all inputs are correct
  ntypes = length(a)
  assert_that(all(c(length(b) == ntypes, length(d) == ntypes, 
                    dim(Xi_as) == c(ntypes, ntypes), dim(Xi_s) == c(ntypes, ntypes),
                    origin_time > 0, origin_type %in% c(0:(ntypes-1)),
                    d >= 0, d < 1, rho > 0, rho <= 1)),
              msg = 'The inputs do not have proper dimensions or values. Please check all parameters!')
  assert_that(all(sapply(seq(1, ntypes), function(i) {all.equal(sum(2*Xi_as[i, ]) + sum(Xi_s[i, ]), 1.)})),
              msg = 'The transition probabilities do not some to 1. Please check!')
  
  # simulate full tree
  tree = sim_adb_origin_complete(origin_time, a, b, d, origin_type, Xi_as, Xi_s, min_tips)
  if (is.null(tree)) {
    return(NULL)
  }
  
  phylo_obj <- prune_tree(tree, rho, min_tips) 
  return(phylo_obj)
}



#' Simulator of the complete Age-Dependent Branching Process for a fixed time interval (since origin)
#' @param origin_time time of birth of the initial particle
#' @param a vector of scale parameters per type
#' @param b vector of shape parameters per type
#' @param d vector of death probabilities per type
#' @param origin_type one of 0,...,n-1 where n is the number of types
#' @param Xi_as matrix of asymetric type transition probabilities
#' @param Xi_s matrix of symetric type transition probabilities
#' @param min_tips minimum number of tips in the phylogeny
#' @export
sim_adb_origin_complete <- function(origin_time, a, b, d, origin_type=0, Xi_as=matrix(0), Xi_s=matrix(1), min_tips = 2) {
  
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
  #nodes = matrix(nrow = 0, ncol = 7, 
  #               dimnames = list(NULL, c("id", "height", "type", "parent", "leftchild", "rightchild", "isTip")))
  
  # sample the lifetime of the first particle
  root_edge = rgamma(1, shape = b[origin_type + 1], scale = a[origin_type + 1])
  nodes = bind_rows(nodes, c(id = 1, height = origin_time - root_edge, type = origin_type, 
                             parent = NA, leftchild = NA, rightchild = NA, status = 1))
  events = nodes
  event_counter = 1
  
  while (nrow(events) > 0) {
    event = as.list(events[1, ]); events = events[-1, ] # look at one event

    if (runif(1) < d[event$type + 1]) {
      # particle dies
      nodes[event$id, "status"] = 0
    } else {
      # particle divides, create two new children
      left_id = event_counter + 1
      right_id = event_counter + 2
      event_counter = event_counter + 2
      nodes[event$id, "status"] = 2

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
      if (event$height - left_lifetime < 0) {
        # censor lifetime
        left_lifetime = event$height
        left_node = c(id = left_id, height = event$height - left_lifetime, type = left_type, 
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      } else {
        left_node = c(id = left_id, height = event$height - left_lifetime, type = left_type, 
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
        events = bind_rows(events, left_node) # to be processed
      }

      right_type = children_types[2]
      right_lifetime = rgamma(1, shape = b[right_type + 1], scale = a[right_type + 1])
      if (event$height - right_lifetime < 0) {
        # censor lifetime, make tip
        right_lifetime = event$height
        right_node = c(id = right_id, height = event$height - right_lifetime, type = right_type, 
                      parent = event$id, leftchild = NA, rightchild = NA, status = 1)
      } else {
        right_node = c(id = right_id, height = event$height - right_lifetime, type = right_type, 
                       parent = event$id, leftchild = NA, rightchild = NA, status = 1)
        events = bind_rows(events, right_node) # to be processed
      }
    
      # add child relationships and append the current edges and edge lengths
      nodes = bind_rows(nodes, left_node, right_node)
      nodes[event$id, c("leftchild", "rightchild")] = c(left_id, right_id)
      edges = rbind(edges, c(event$id, left_id), c(event$id, right_id))
      edge_lengths = c(edge_lengths, left_lifetime, right_lifetime)
    }
  }
  
  # check number of tips
  Nnode = sum(nodes$status == 2)
  Ntip = sum(nodes$status %in% c(0,1))
  if (Ntip < min_tips) {
    message("The simulated tree has too few tips. Try another seed.")
    return(NULL)
  }
  
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
  tree@phylo$origin = origin_time
  data = as_tibble(tree)
  types = nodes %>% 
    select(node = label, status, type) %>% 
    mutate(type = as.factor(type)) %>%
    arrange(node) %>%
    as_tibble
  tree@data = types
  
  return(tree)
}

