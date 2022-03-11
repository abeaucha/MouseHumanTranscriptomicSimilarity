# HumanCorticalLabels.R -------------------------------------------------------
#
# Script to extract human cortical labels from expression matrix 
# and save as a separate file
#
# Antoine Beauchamp
# Created: February 7th, 2022
# Edited: February 7th, 2022

#Libraries
library(tidyverse)
library(data.tree)

#Tree functions
source("../../../functions/tree_tools.R")

#Import human sample expression matrix for labels
dfExprHuman <- suppressMessages(read_csv("../../../data/HumanExpressionMatrix_Samples_pipeline_v1_labelled.csv"))

#Import human label sets
load("../../../data/TreeLabels.RData")

#Import human tree
load("../../../AHBA/data/HumanExpressionTree.RData")
treeHuman <- Clone(treeHumanExpr)
rm(treeHumanExpr)

#Get all cortical leaf nodes
corticalRegionsAll <- FindNode(treeHuman, "cerebral cortex")$Get("name", filterFun = isLeaf)

#Create a copy of the tree to prune down to be bilateral
treeHumanBilateral <- Clone(treeHuman)

#Remove white matter and ventricles from the tree
pruneAnatTree(treeHumanBilateral, nodes = c("white matter", "sulci & spaces"), method = "AtNode")

#Prune the tree down to 166 bilateral regions
pruneAnatTree(treeHumanBilateral, nodes = listLabelsHuman$Region166, method = "BelowNode")

#Extract cortical regions from the tree
corticalRegionsBilateral <- FindNode(treeHumanBilateral, "cerebral cortex")$Get("name", filterFun = isLeaf)

#Empty data frame to match lateral leaf nodes with bilateral parent regions
dfCorticalRegionsAll <- tibble(RegionAll = corticalRegionsAll,
                               Region166 = "")

#Match leaf nodes with bilateral parents
for (i in 1:nrow(dfCorticalRegionsAll)){
  cxRegion <- dfCorticalRegionsAll$RegionAll[[i]]
  cxRegionPath <- FindNode(treeHuman, cxRegion)$path
  indBilateral <- which(cxRegionPath %in% corticalRegionsBilateral)
  dfCorticalRegionsAll[i,"Region166"] <- cxRegionPath[indBilateral]  
}

#Combine cortical labels at multiple levels with lateral leaf nodes
dfCorticalRegions <- dfExprHuman %>%
  select(Region16, Region88, Region166) %>% 
  distinct() %>% 
  filter(Region166 %in% corticalRegionsBilateral) %>% 
  left_join(dfCorticalRegionsAll, by = "Region166") %>% 
  mutate(RegionAll = factor(RegionAll, levels = corticalRegionsAll)) %>% 
  arrange(RegionAll)

#Save cortical labels as csv file
write_csv(x = dfCorticalRegions,
          file = "HumanCorticalLabels.csv")
   