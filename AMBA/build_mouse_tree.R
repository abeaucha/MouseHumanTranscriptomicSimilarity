# ----------------------------------------------------------------------------
# build_mouse_tree.R
# Author: Antoine Beauchamp
# Created: January 31st, 2022
#
# Build a gene expression tree
#
# Description
# -----------
# This script maps the values from a voxel-wise expression matrix to the
# Allen Mouse Brain Atlas hierarchical ontology. 


# Packages -------------------------------------------------------------------

message("Initializing...")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(rjson))


# Functions ------------------------------------------------------------------

#Import tree functions
source("../functions/tree_tools.R")


# Variables ------------------------------------------------------------------

dirData <- "data/"
dirImaging <- "data/imaging/"
fileExpr <- "MouseExpressionMatrix_ROI_DSURQE_coronal_maskcoronal_log2_grouped_imputed.csv"
fileTree <- "DSURQE_tree.json"
fileAtlasDefs <- "DSURQE_40micron_R_mapping_long.csv"


# Import ---------------------------------------------------------------------

#Import expression data
message("Importing expression data...")
dfExpr <- suppressMessages(read_csv(str_c(dirData, fileExpr)))

message("Generating expression tree...")

#Import DSURQE/AMBA tree from JSON
treeDefs <- parse_abi_hierarchy(str_c(dirData, fileTree))
treeMouseExpr <- Clone(treeDefs)

#Import DSURQE atlas definitions
defs <- suppressMessages(read_csv(str_c(dirImaging, fileAtlasDefs)))


# Assign values to tree ------------------------------------------------------

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


# Write ----------------------------------------------------------------------

message("Writing to file...")

fileOut <- "MouseExpressionTree_DSURQE.RData"

save(treeMouseExpr, 
     file = str_c(dirData, fileOut))
