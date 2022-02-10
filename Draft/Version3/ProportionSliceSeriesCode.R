library(tidyverse)
library(data.tree)
library(RMINC)
library(MRIcrotome)

source("/projects/abeauchamp/Functions/TreeTools.R")

#Load tree labels and tree
load("/projects/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/Data/TreeLabelsReordered.RData")
load("/projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/MouseExpressionTree_DSURQE.RData")
treeMouse <- Clone(treeMouseExpr)
rm(treeMouseExpr)

#Remove white matter and ventricles
pruneAnatTree(treeMouse,
              nodes = c("fiber tracts", "ventricular systems"),
              method = "AtNode")

#Import DSURQE images in AMBA space, 200um
dsurqeLabels_200um <- mincGetVolume("/projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/MRI/DSURQE_Allen_labels.mnc")
dsurqeMask_200um <- mincGetVolume("/projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/MRI/coronal_200um_coverage_bin0.8.mnc")
dsurqeAverage_200um <- mincGetVolume("/projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/MRI/DSURQE_Allen_average.mnc")

#Create a correspondence between voxels in the image space and in the expression matrices

#Indices for the mask and for grey matter ROIs
indMask <- dsurqeMask_200um == 1
indGM <- dsurqeLabels_200um %in% treeMouse$`Basic cell groups and regions`$label
indVoxels <- indMask & indGM

#Create an empty array and assign IDs to non-null voxels
emptyArray <- numeric(length(dsurqeLabels_200um))
names(emptyArray) <- character(length(emptyArray))
voxelNames <- str_c("V", 1:sum(indVoxels))
names(emptyArray)[indVoxels] <- voxelNames

#Data frame containing proportions for striatal voxels
load("/projects/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/Draft/Version3/StriatumProportions.RData")

#Initialize 200um MINC arrays
listArrayHumanTargets_200um <- vector(mode = "list", length = ncol(dfSimStriatum_MaxSim_Proportions_HumanTargets))
names(listArrayHumanTargets_200um) <- colnames(dfSimStriatum_MaxSim_Proportions_HumanTargets)
listArrayHumanTargets_200um <- map(listArrayHumanTargets_200um, function(x){return(emptyArray)})

#For each of the human targets, map the proportions to voxels in a MINC array
for(j in 1:ncol(dfSimStriatum_MaxSim_Proportions_HumanTargets)){
  indVoxelMatch <- match(rownames(dfSimStriatum_MaxSim_Proportions_HumanTargets), names(emptyArray))
  listArrayHumanTargets_200um[[j]][indVoxelMatch] <- dfSimStriatum_MaxSim_Proportions_HumanTargets[,j]
  attributes(listArrayHumanTargets_200um[[j]]) <- attributes(dsurqeLabels_200um)
}

#Files and paths for resampling
outfiles <- str_c("Figure6_ss_", str_remove(names(listArrayHumanTargets_200um), " "), "_200um.mnc")
outpaths <- str_c("/projects/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/Draft/Version3/", outfiles)
anatomy50um <- "/projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/MRI/average_template_50um.mnc"
newpaths <- str_replace(outpaths, "200um", "50um")
resampleCommands <- str_c("mincresample",
                          "-like",
                          anatomy50um,
                          outpaths,
                          newpaths,
                          "-clobber",
                          sep = " ")

#Initialize 50um arrays
listArrayHumanTargets_50um <- vector(mode = "list", length(listArrayHumanTargets_200um))
names(listArrayHumanTargets_50um) <- names(listArrayHumanTargets_200um)

#Resample MINC arrays and reload
for(i in 1:length(listArrayHumanTargets_200um)){
  
  mincWriteVolume(listArrayHumanTargets_200um[[i]],
                  output.filename = outpaths[i],
                  like.filename = attributes(listArrayHumanTargets_200um[[i]])$likeVolume,
                  clobber = TRUE)
  
  system(resampleCommands[i])
  
  listArrayHumanTargets_50um[[i]] <- mincGetVolume(newpaths[i])
  
}

dsurqeMask_50um <- mincGetVolume("/projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/MRI/average_template_50um_mask.mnc")
dsurqeAverage_50um <- mincGetVolume("/projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/MRI/DSURQE_Allen_average_50um.mnc")
dsurqeAverage_50um[dsurqeMask_50um == 0] <- 0

#Overlay range
proplow = 0.1
prophigh = 1

#200um slice series
sliceSeries(nrow = 1, ncol = 8, begin = 36, end = 49) %>% 
  anatomy(mincArray(dsurqeAverage_200um), low = 700, high = 1400) %>%
  overlay(mincArray(listArrayHumanTargets_200um[["caudate nucleus"]]), low = proplow, high = prophigh) %>%
  addtitle("caudate") %>% 
  sliceSeries() %>% anatomy() %>%
  overlay(mincArray(listArrayHumanTargets_200um[["putamen"]]), low = proplow, high = prophigh) %>%
  addtitle("putamen") %>% 
  sliceSeries() %>% anatomy() %>%
  overlay(mincArray(listArrayHumanTargets_200um[["nucleus accumbens"]]), low = proplow, high = prophigh) %>%
  addtitle("accumbens") %>% 
  sliceSeries() %>% anatomy() %>%
  overlay(mincArray(listArrayHumanTargets_200um[["septal nuclei"]]), low = proplow, high = prophigh) %>%
  addtitle("septal nuclei") %>% 
  draw(layout = "row")

#50um slice series
sliceSeries(nrow = 1, ncol = 8, begin = 144, end = 200) %>%
  anatomy(mincArray(dsurqeAverage_50um), low = 700, high = 1400) %>%
  overlay(mincArray(listArrayHumanTargets_50um[["caudate nucleus"]]), low = proplow, high = prophigh) %>%
  sliceSeries() %>% anatomy() %>%
  overlay(mincArray(listArrayHumanTargets_50um[["putamen"]]), low = proplow, high = prophigh) %>%
  sliceSeries() %>% anatomy() %>%
  overlay(mincArray(listArrayHumanTargets_50um[["nucleus accumbens"]]), low = proplow, high = prophigh) %>%
  sliceSeries() %>% anatomy() %>%
  overlay(mincArray(listArrayHumanTargets_50um[["septal nuclei"]]), low = proplow, high = prophigh) %>% 
  draw(layout = "row")
