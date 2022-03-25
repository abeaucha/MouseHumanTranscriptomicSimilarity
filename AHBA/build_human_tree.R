# ----------------------------------------------------------------------------
# build_human_tree.R 
# Author: Antoine Beauchamp
# Created: March 9th, 2021
# 
# Build the gene expression tree for the Allen Human Brain Atlas 


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
  make_option("--infile",
              type = "character",
              help = paste("Name of CSV file containing the sample-wise",
                           "expression matrix.")),
  make_option("--outfile",
              type = "character",
              default = "HumanExpressionTree.RData",
              help = paste("Name of .RData file in which to export the tree.",
                           "[default %default]")),
  make_option("--samplefile",
              type = "character",
              help = "Name of CSV file containing sample information."),
  make_option("--treefile",
              type = "character",
              help = paste("Name of JSON file containing hierarchical ontology.",
                           "[default %default]"))
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args[["infile"]])){
  stop("Argument --infile empty with no default.")
}

if (is.null(args[["samplefile"]])){
  stop("Argument --samplefile empty with no default.")
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


# Map microarray samples to tree nodes ---------------------------------------

dirData <- args[["datadir"]]
fileExpr <- args[["infile"]]
fileSampleInfo <- args[["samplefile"]]
fileTreeDefs <- args[["treefile"]]

message("Mapping microarray samples to the AHBA ontology...")

#Load AHBA sample information
dfSampleInfo <- suppressMessages(read_csv(file.path(dirData,fileSampleInfo)))

#Load the human tree definitions
treeHumanDefs <- parse_abi_hierarchy(file.path(dirData, fileTreeDefs))

#Make cerebellar vermis names different from cerebellar hemispheres
#Otherwise this messes with the sample assignment
treeHumanDefs$Do(function(node){
  if(isLeaf(node)){
    if("vermis" %in% node$path){
      node$name <- str_c("vermal ", node$name)
    }
  }
})


dfSampleInfo <- dfSampleInfo %>% 
  mutate(structure_name = ifelse(str_detect(structure_acronym, "^Ve"), 
                                 str_c("vermal ", structure_name), 
                                 structure_name)) 

#Assign AHBA samples to nodes on the tree
treeHumanDefs$Do(function(node){
  if(node$name %in% dfSampleInfo$structure_name){
    whichInd <- which(node$name == dfSampleInfo$structure_name)
    node$samples <- dfSampleInfo$SampleID[whichInd]
  }
})

#Aggregate samples up the tree
treeHumanDefs$Do(function(node){
  node$samples <- c(node$samples, unlist(Aggregate(node, "samples", c)))
}, traversal = "post-order")

#The aggregation generates duplicate samples and NAs. Remove those
treeHumanDefs$Do(function(node){
  node$samples <- unique(node$samples[!is.na(node$samples)])
})

#Remove nodes that don't have any samples
Prune(treeHumanDefs, pruneFun = function(node){length(node$samples) != 0})


# Add gene expression values to the tree -------------------------------------

message("Adding gene expression data to the tree...")

#Import sample expression matrix
matHumanExpr <- suppressMessages(read_csv(str_c(file.path(dirData, fileExpr)))) %>% 
  column_to_rownames("Gene") %>% 
  as.matrix()

#Aggregate sample expression values at each node
treeHumanExpr <- Clone(treeHumanDefs)
treeHumanExpr$Do(function(node){
  
  sampleInd <- match(node$samples, colnames(matHumanExpr))
  
  if(length(sampleInd) != 0){
    node$Expression <- rowMeans(as.matrix(matHumanExpr[,sampleInd]))
  } else {
    node$Expression <- rep(0, nrow(matHumanExpr))
  }
  
  node$Gene <- rownames(matHumanExpr)
  
})


# Write to file --------------------------------------------------------------

message("Writing data to file...")

save(treeHumanExpr,
     file = file.path(dirData, args[["outfile"]]))

