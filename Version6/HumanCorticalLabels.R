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
source("/projects/abeauchamp/Functions/TreeTools.R")

#Import human sample expression matrix for labels
dfExprHuman <- suppressMessages(read_csv("/projects/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/Data/HumanExpressionMatrix_Samples_labelled.csv"))

#Import human label sets
load("/projects/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/Data/TreeLabels.RData")

#Import human tree
load("/projects/abeauchamp/Projects/MouseHumanMapping/AllenHumanBrainAtlas/Data/HumanExpressionTree_Samples.RData")
treeHuman <- Clone(treeHumanExpr)
rm(treeHumanExpr)

#Remove white matter and ventricles from the tree
pruneAnatTree(treeHuman, nodes = c("white matter", "sulci & spaces"), method = "AtNode")

#Prune the tree down to 166 bilateral regions
pruneAnatTree(treeHuman, nodes = listLabelsHuman$Region166, method = "BelowNode")

#Extract cortical regions from the tree
corticalRegions <- FindNode(treeHuman, "cerebral cortex")$Get("name", filterFun = isLeaf)

#Filter sample labels based on cortical regions
dfCorticalRegions <- dfExprHuman %>% 
  select(Region16, Region88, Region166) %>% 
  distinct() %>% 
  filter(Region166 %in% corticalRegions) %>% 
  mutate(Region166 = factor(Region166, levels = corticalRegions)) %>% 
  arrange(Region166)

#Save cortical labels as csv file
write_csv(x = dfCorticalRegions,
          file = "/projects/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/Draft/Version6/HumanCorticalLabels.csv")
   