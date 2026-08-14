#' Ornstein-Uhlenbeck Model for Gene Expression along tree
#' @param tree a complete tree (with discrete types)
#' @param n_trait total number of traits to sim
#' @param theta matrix where mu_ij is the mean of trait i in type j
#' @param sigma matrix where sigma_ij is the sd of trait i in type j
#' @param alpha matrix where alpha_ij is the strength of selection on trait i in type j

sim_expression_ontree <- function(tree,n_trait,n_type,theta,sigma,alpha){
  
  
  return(NULL)
  
}

sim_expression_ontree_gene <- function(tree, 
                            alpha = 0.5, 
                            sigmasq = 0.2,
                            root_value = 0,
                            theta = NULL # df with columns type and theta 
){
  
  
  sigma = sqrt(sigmasq)
  
  evolved_expression = tree %>% as_tibble() %>% as.data.frame()
  root = evolved_expression %>% filter(parent == node) %>% .$node
  
  if(nrow(theta) != length(unique(tree@data$type))){
    message("Error: the number of thetas provided does not match the number of types in the tree. Provide the thetas as a data frame with columns 'type' and 'theta'")
    return(NULL)
  }
  
  evolved_expression$expr = NA
  evolved_expression$expr[which(evolved_expression$node== root)] = root_value
  
  evolved_expression_noise = evolved_expression %>%
    filter(node != root) %>% 
    rowwise() %>%
    mutate(noise = case_when(
      alpha == 0 ~ rnorm(1, # reduces to BM
                         mean = 0, 
                         sd = sigma * sqrt(branch.length)),
      TRUE ~ rnorm(1, 
                   mean = 0,
                   sd = sigma * sqrt((1 - exp(-2 * alpha * branch.length)) / (2 * alpha)))
    )) %>%
    ungroup() %>%
    select(node, noise)
  
  evolved_expression = evolved_expression %>% 
    left_join(evolved_expression_noise, by = c("node")) %>%# if I apply directly rowwise to tree_tibble, I lose the treedata abstraction
    left_join(thetas, by = c("type"))
  
  evolved_expression$noise[which(evolved_expression$node== root)] = 0 
  
  
  while(anyNA(evolved_expression$expr)){
    
    evolved_expression = evolved_expression %>%
      dplyr::select(-c(any_of(c("parent_expr")))) %>% # to allow overwriting
      left_join(evolved_expression %>% 
                  dplyr::select(parent = node,
                                parent_expr = expr),
                by = "parent") %>%
      mutate(expr = case_when(!is.na(expr) ~ expr,
                              is.na(parent_expr) ~ expr,
                              T ~ theta + (parent_expr-theta)*exp(-alpha*branch.length) + noise 
      )) 
    
  }
  
  
  evolved_expression %<>% select(node, type, expr)
  
  return(evolved_expression)
  
}

