# build_mouse_tree.R ------------------------------------------------------------------
# 
#
# Antoine Beauchamp
# Created: January 31st, 2022
# Edited March 8th, 2022
# ---------------------------------------------------------------------------------------


message("Initializing...")

#Libraries
library(tidyverse)
library(data.tree)
library(rjson)


#Import tree functions
source("../functions/abi_parsing_funs.R")

#Import expression data

message("Importing expression data...")

# dfExpr <- read_csv("data/<file>")
dfExpr <- read_csv(str_c(pathAMBA, "Data/", "MouseExpressionMatrix_ROI_DSURQE_coronal_maskcoronal.csv"))

message("Generating expression tree...")

#Import DSURQE/AMBA tree from JSON
treeDefs <- parse_abi_hierarchy("data/DSURQE_tree.json")
treeMouseExpr <- Clone(treeDefs)

#Import DSURQE atlas definitions
defs <- read_csv("data/imaging/DSURQE_40micron_R_mapping_long.csv")

#Assign DSURQE ROI expression values to leaf nodes on the tree
treeMouseExpr$Do(function(node){
  if(isLeaf(node)){
    whichInd <- which(defs$Label == node$label)
    node$Expression <- dfExpr[,colnames(dfExpr) == defs$Structure[whichInd]][[1]]
  }
})

#Aggregate expression values up the tree
treeMouseExpr$Do(function(node){
  node$Expression <- Aggregate(node, "Expression", rowMeans)
}, traversal = "post-order")


#Assign corresponding gene names to the tree
treeMouseExpr$Do(function(node){
  node$Gene <- dfExpr$Gene
})


message("Writing to file...")

save(treeMouseExpr, 
     file = "data/MouseExpressionTree_DSURQE.RData"))

