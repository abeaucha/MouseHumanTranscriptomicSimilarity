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


# Functions ------------------------------------------------------------------

source("../functions/tree_tools.R")


# Map microarray samples to tree nodes ---------------------------------------

message("Mapping microarray samples to the AHBA ontology...")

#Load AHBA sample information
fileSampleInfo <- "SampleInformation_pipeline_v1.csv"
dfSampleInfo <- suppressMessages(read_csv(paste0("data/",fileSampleInfo)))

#Load the human tree definitions
fileTreeDefs <- "AHBA_hierarchy_definitions.json"
treeHumanDefs <- parse_abi_hierarchy(paste0("data/", fileTreeDefs))

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
fileExpr <- "HumanExpressionMatrix_Samples_pipeline_v1.csv"
matHumanExpr <- suppressMessages(read_csv(str_c(paste0("data/", fileExpr)))) %>% 
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

fileOut <- "HumanExpressionTree.RData"
save(treeHumanExpr,
     file = paste0("data/", fileOut))

