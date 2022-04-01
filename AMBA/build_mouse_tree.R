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

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(optparse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--datadir",
              type = "character",
              default = "data/",
              help = paste("Directory containing expression matrix CSV files.",
                           "[default %default]")),
  make_option("--imgdir",
              type = "character",
              default = "data/imaging/",
              help = "Directory containing imaging files. [default %default]"),
  make_option("--infile",
              type = "character",
              help = paste("Name of CSV file containing the region-wise",
                           "expression matrix.")),
  make_option("--outfile",
              type = "character",
              default = "MouseExpressionTree_DSURQE.RData",
              help = paste("Name of .RData file in which to export the tree.",
                           "[default %default]")),
  make_option("--treefile",
              type = "character",
              default = "DSURQE_tree.json",
              help = paste("Name of JSON file containing hierarchical", 
                           "ontology. [default %default]")),
  make_option("--defs",
              type = "character",
              help = paste("Name of CSV file containing the names of the", 
                           "neuroanatomical regions corresponding to the",
                           "nodes of the tree in --treefile")),
  make_option("--verbose",
              type = 'character',
              default = 'true',
              help = "Verbosity. [default %default]")
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args[["infile"]])){
  stop("Argument --infile empty with no default.")
}

if (is.null(args[["defs"]])){
  stop("Argument --defs empty with no default.")
}


# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset("--file=") %>% 
  str_remove("--file=") %>%
  dirname()  

path_tree_tools <- str_c(working_dir, 
                         script_dir, 
                         "../functions/tree_tools.R", 
                         sep = "/")

source(path_tree_tools)


# Import ---------------------------------------------------------------------

dirData <- args[["datadir"]]
dirImaging <- args[["imgdir"]]
fileExpr <- args[["infile"]]
fileTree <- args[["treefile"]]
fileAtlasDefs <- args[["defs"]]
verbose <- ifelse(args[['verbose']] == 'true', TRUE, FALSE)

#Import expression data
if (verbose) {message("Importing expression data...")}
dfExpr <- suppressMessages(data.table::fread(str_c(dirData, fileExpr),
                                             header = TRUE)) %>% 
  as_tibble()

if (verbose) {message("Generating expression tree...")}

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

if (verbose) {message(str_c("Writing to file: ", args[["outfile"]], "..."))}

save(treeMouseExpr, 
     file = str_c(dirData, args[["outfile"]]))
