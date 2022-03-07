# BuildHumanExprTree --------------------------------------------------------------------
# 
#
#
# Antoine Beauchamp
# Created: March 9th, 2021
# Edited: August 23rd, 2021
# ---------------------------------------------------------------------------------------

# Libraries -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
library(data.tree)
library(rjson)


# Functions -----------------------------------------------------------------------------

source("/projects/abeauchamp/Functions/abi_parsing_funs.R")
source("/projects/abeauchamp/Functions/TreeTools.R")


# Map microarray samples to tree nodes --------------------------------------------------

message("Mapping AHBA microarray samples to the human neuroanatomical tree...")

#Paths
pathAHBA <- "/projects/abeauchamp/Projects/MouseHumanMapping/AllenHumanBrainAtlas/"

#Load AHBA sample information
dfSampleInfo <- suppressMessages(read_csv(str_c(pathAHBA, "Data/", "SampleInformation.csv")))

#Load the human tree definitions
treeHumanDefs <- parse_abi_hierarchy("/projects/abeauchamp/Atlases/AHBA/human_hierarchy_definitions.json")

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


# Add gene expression values to the tree ------------------------------------------------

message("Adding gene expression data to the tree...")

#Import sample expression matrix
matHumanExpr <- suppressMessages(read_csv(str_c(pathAHBA, "Data/", "HumanExpressionMatrix_Samples.csv"))) %>% 
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


# Write to file -------------------------------------------------------------------------

message("Writing data to file...")

save(treeHumanExpr,
     file = str_c(pathAHBA, "Data/", "HumanExpressionTree_Samples.RData"))

