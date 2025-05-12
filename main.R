#####################################
# --MAIN MODEL DEFINITION & SETUP---#
#####################################
# Simulator(s) based on cell lineage trees
# ------------>
#         ┌── o
#     ┌───┤
#     │   └── o
# ───┤
#     │   ┌── o
#     └───┤
#         └── o
#
# load packages
#source("~/ADBP/R/simulate_adbp_ntaxa.R") #as it's own package eventually..? or incorporate?
source("~/ADBP/R/simulate_adbp_origin.R")
source("models/gene_expression.R")
library(ape)


simulate_cellular_processes <- function(
    fname, #file name for output
    origin_time,#origin time of process
    a,#scale
    b=1,#shape
    d=0,#death prob
    samp=1, #sampling prob (at present only)
    origin_type=0,#cell type at origin
    Xi_as=matrix(0),#asymmetric type changes
    Xi_s=matrix(0),#symmetric type changes
    simulate_expression=F,
    lineage_tracing=F,
    output_dir = "./output"
    ){

    # LINEAGE TREE (COMPLETE TREE)
    complete_tree <- simulate_complete_tree(
        origin_time,a,b,d,origin_type,
        Xi_as, Xi_s)

    # GENE EXPRESSION DATA
    if (simulate_expression) {
        # model returns 
        G <- developmental_OU_model(
            complete_tree)
        # add error model?
    }
    
    # LINEAGE TRACING DATA
    if (lineage_tracing){
        #--to do
    }
    
    
    # SAMPLE TIPS:
    
    # -----OUTPUT-----#
    # 1. output lineage tree (newick)

    # 2. data
}
        
