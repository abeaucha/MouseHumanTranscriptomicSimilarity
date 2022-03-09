# ProcessLabelledData.R
#
# Antoine Beauchamp
# Created: August 25th, 2021
# Edited: August 25th, 2021

#Libraries
suppressPackageStartupMessages(library(tidyverse))

#Functions
source("/projects/abeauchamp/Projects/MouseHumanMapping/Functions/ProcessingTools.R")

#Set path
pathHome <- "/projects/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/"

#Import data
dfExprMouse <- suppressMessages(read_csv(str_c(pathHome, "Data/", "MouseExpressionMatrix_Voxel_coronal_maskcoronal_imputed_labelled.csv")))
dfExprHuman <- suppressMessages(read_csv(str_c(pathHome, "Data/", "HumanExpressionMatrix_Samples_labelled.csv")))

#Extract labels from data frames
dfLabelsMouse <- dfExprMouse %>% select(contains("Region"))
dfLabelsHuman <- dfExprHuman %>% select(contains("Region"))

#Normalize mouse data
dfExprMouse_scaled <- dfExprMouse %>% 
  select(-contains("Region")) %>% 
  as.matrix() %>% 
  scaler(axis = "rows") %>% 
  scaler(scale = FALSE, axis = "columns") %>% 
  as_tibble() %>% 
  bind_cols(dfLabelsMouse)

#Normalize human data
dfExprHuman_scaled <- dfExprHuman %>% 
  select(-contains("Region")) %>% 
  as.matrix() %>% 
  scaler(axis = "rows") %>% 
  scaler(scale = FALSE, axis = "columns") %>% 
  as_tibble() %>% 
  bind_cols(dfLabelsHuman)

#Extract genes list from mouse data (same as human)
genes <- colnames(dfExprMouse_scaled)[!str_detect(colnames(dfExprMouse_scaled), "Region")]


mouseLabels <- "Region67"
humanLabels <- "Region88"

#Aggregate mouse expression data under label set
dfExprMouse_scaled_means <- dfExprMouse_scaled %>% 
  select(Region = mouseLabels, genes) %>% 
  group_by(Region) %>% 
  summarise_all(mean) %>% 
  ungroup()

#Aggregate human expression data under label set
dfExprHuman_scaled_means <- dfExprHuman_scaled %>% 
  select(Region = humanLabels, genes) %>% 
  group_by(Region) %>% 
  summarise_all(mean) %>% 
  ungroup()

#Write scaled mouse voxel data to file
write_csv(dfExprMouse_scaled,
          path = str_c(pathHome, "Data/", "MouseExpressionMatrix_Voxel_coronal_maskcoronal_imputed_labelled_scaled.csv"))

#Write scaled human sample data to file
write_csv(dfExprHuman_scaled,
          path = str_c(pathHome, "Data/", "HumanExpressionMatrix_Samples_labelled_scaled.csv"))
  
#Write aggregated mouse data to file
write_csv(dfExprMouse_scaled_means,
          path = str_c(pathHome, "Data/", "MouseExpressionMatrix_ROI_", mouseLabels, ".csv"))

#Write aggregated human data to file
write_csv(dfExprHuman_scaled_means,
          path = str_c(pathHome, "Data/", "HumanExpressionMatrix_ROI_", humanLabels, ".csv"))
