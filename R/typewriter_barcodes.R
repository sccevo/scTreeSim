
colouring <- function(x,l, lambda,sample_p,chars){
	p1 <- stringr::str_split(x,"")[[1]]
    t <- 0
    M <- length(p1)
    while (t <= l){
        time_to_mut <- rexp(1,lambda)
        t <- t + time_to_mut	
        if (p1[M] == "0" && t <= l) {
            pos <- min(which(p1 == "0"))
            new <- sample(chars, size=1, prob=sample_p) 
            p1[pos] <- new
        } 
    }
    return(paste(p1, collapse=""))
}


sim_typewriter_barcode <- function(rate, chars, sampling, tree,root){
	# first, mutate along root edge	
	new_root <- colouring(root, tree$root.edge, rate, sampling, chars)	
	mut_tree <- rTraitMult(tree, 
        	colouring, 
        	root.value=new_root, 
        	ancestor=TRUE, 
        	lambda=rate, 
        	sample_p=sampling, 
        	chars=chars)
   	ntip <- tree$Nnode+1
  	Y <- mut_tree$x1[1:ntip]
    return(Y)
}

#' Simulator of Sequentially-Edited Barcodes
#' @param tree a time-scaled phylogeny
#' @param k number of target sites
#' @param lambda_vec array of barcode-specific editing rates (length m)
#' @param m number of barcodes per cell (int >= 1)
#' @param chars array of unique characters to be insterted
#' @export
typewriter_barcodes <- function(tree, k, lambda_vec, m, chars){
	# sim bars over tree
	root <- paste0(rep("0",k), collapse='')
    sample_p <- rep(1/length(chars), length(chars))
	res <- lapply(lambda_vec, sim_typewriter_barcode, chars, sample_p,tree,root)
	return(res)
}
