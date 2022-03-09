# LabelVoxelData.R ------------------------------------------------------------
#
# 
# 
#
# Antoine Beauchamp
# Edited: August 23rd, 2021

# Libraries -------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
library(RMINC)
library(data.tree)
library(optparse)


# Command line arguments ------------------------------------------------------

option_list <- list(
  make_option("--mousedata",
              type = "character",
              default = "MouseExpressionMatrix_Voxel_coronal_maskcoronal_imputed.csv",
              help = "Mouse CSV file containing voxel data to label. File must exist at /projects/abeauchamp/Projects/MouseHumanMapping/AllenMouseBrainAtlas/Data/ [default %default]"),
  make_option("--humandata",
              type = "character",
              default = "HumanExpressionMatrix_Samples.csv",
              help = "Human CSV file containing sample data to label. File must exist at /projects/abeauchamp/Projects/MouseHumanMapping/AllenHumanBrainAtlas/Data/ [default %default]"),
  make_option("--savemouse",
              type = "character",
              default = "true",
              help = "Option to save labelled mouse data to file [default %default]"),
  make_option("--savehuman",
              type = "character",
              default = "true",
              help = "Option to save labelled human data to file [default %default]")
)

args = parse_args(OptionParser(option_list = option_list))

if (!(args[["savemouse"]] %in% c("true", "false"))) {
  stop(str_c("Argument --savemouse must be one of [true, false] (got ", args[["savemouse"]], ")"))
}

if (!(args[["savehuman"]] %in% c("true", "false"))) {
  stop(str_c("Argument --savehuman must be one of [true, false] (got ", args[["savehuman"]], ")"))
}



# Functions -------------------------------------------------------------------

source("/projects/abeauchamp/Functions/TreeTools.R")
source("/projects/abeauchamp/Projects/MouseHumanMapping/Functions/ProcessingTools.R")

#' Label voxels/samples with neuroanatomical regions
#'
#' @param measurements 
#' @param tree 
#' @param treefield 
#'
#' @return
labelRegions <- function(measurements, tree, treefield){
  
  #Get a list of all measurements and their corresponding ROIs
  listMeasurements <- tree$Get(treefield, filterFun = isLeaf)
  vecMeasurements <- unlist(listMeasurements)
  
  #Unlisting to vector will break the names. Fix them
  structNames <- names(listMeasurements)
  measurementNames <- names(vecMeasurements)
  
  for (struct in structNames){
    
    structRgx <- struct %>% 
      str_replace("\\(", "\\\\(") %>% 
      str_replace("\\)", "\\\\)")
    
    structRgx <- str_c("^", structRgx, "[0-9]+", "$")
    
    indStruct <- str_which(measurementNames, structRgx)
    
    measurementNames[indStruct] <- struct
  }
  
  #Proper names  
  names(vecMeasurements) <- measurementNames

  indMatchStructs <- match(measurements, vecMeasurements)
  
  measurementStructs <- names(vecMeasurements)[indMatchStructs]
  
  return(measurementStructs)
}


# Paths -----------------------------------------------------------------------

#Set paths
pathHome <- "/projects/abeauchamp/Projects/MouseHumanMapping/"
pathMouse <- str_c(pathHome, "AllenMouseBrainAtlas/Data/")
pathHuman <- str_c(pathHome, "AllenHumanBrainAtlas/Data/")

setwd(str_c(pathHome, "Paper_Descriptive"))

fileMouse <- args[["mousedata"]]
fileHuman <- args[["humandata"]]

if (!(fileMouse %in% list.files(pathMouse))) {
  stop("Mouse file ", fileMouse, " not found at ", pathMouse)
} else if (!(fileMouse %in% str_subset(list.files(pathMouse), "^MouseExpressionMatrix_Voxel"))) {
  stop("Mouse file ", fileMouse, " does not contain voxelwise expression data")
}

if(!(fileHuman %in% list.files(pathHuman))){
  stop("Human file ", fileHuman, " not found at ", pathHuman)
} else if (!(fileHuman %in% str_subset(list.files(pathHuman), "^HumanExpressionMatrix_Samples"))) {
  stop("Human file ", fileHuman, " does not contain samplewise expression data")
}

message("Labelling data from mouse file: ", fileMouse)
message("Labelling data from human file: ", fileHuman)

maskFlag <- str_extract(fileMouse, "mask[a-z]+")

pathFileMouse <- str_c(pathMouse, fileMouse)
pathFileHuman <- str_c(pathHuman, fileHuman)


# Importing and processing ----------------------------------------------------

message("Importing data...")

#Load expression data
dfExprMouse <- suppressMessages(read_csv(pathFileMouse))
dfExprHuman <- suppressMessages(read_csv(pathFileHuman))

#Subset genes for mouse-human homologs
listExpr <- intersectGeneHomologs(data = list(Mouse = dfExprMouse, 
                                              Human = dfExprHuman),
                                  homologs = str_c(pathHome, "Paper_Descriptive/Data/MouseHumanGeneHomologs.csv"))

#Extract the data frames from list
dfExprMouse <- listExpr$Mouse
dfExprHuman <- listExpr$Human
rm(listExpr)

#Extract genes and remove from df
genesMouse <- dfExprMouse$Gene
genesHuman <- dfExprHuman$Gene

dfExprMouse <- dfExprMouse %>% select(-Gene)
dfExprHuman <- dfExprHuman %>% select(-Gene)

#Clean up mouse column names
colnames(dfExprMouse) <- str_c("V", colnames(dfExprMouse))



# Build the mouse data tree ---------------------------------------------------

message("Importing mouse data tree...")

#Load DSURQE atlas and mask
dsurqe <- mincGetVolume(str_c(pathMouse, "MRI/DSURQE_Allen_labels.mnc"))

if (maskFlag == "masksagittal"){
  mask <- mincGetVolume(str_c(pathMouse, 'MRI/sagittal_200um_coverage_bin0.8.mnc'))
} else if (maskFlag == "maskcoronal") {
  mask <- mincGetVolume(str_c(pathMouse, 'MRI/coronal_200um_coverage_bin0.8.mnc'))
} else {
  stop(str_c("Invalid maskFlag value: ", maskFlag))
}

#Mask the DSURQE atlas
dsurqe_masked <- dsurqe[mask == 1]

#Load the mouse tree
load(str_c(pathMouse, 'MouseExpressionTree_DSURQE.RData'))

treeMouse <- Clone(treeMouseExpr)
rm(treeMouseExpr)

#Assign voxels to leaf nodes on the tree
treeMouse$Do(function(node){
  if(isLeaf(node)){
    node$voxels <- colnames(dfExprMouse)[dsurqe_masked == node$label]
  }
})

#Aggregate voxel names up the tree
treeMouse$Do(function(node){
  node$voxels <- unlist(Aggregate(node, "voxels", c))
})

#Remove white matter and ventricles
cutAtNodes <- c("fiber tracts", "ventricular systems")
pruneAnatTree(treeMouse, nodes = cutAtNodes, method = "AtNode")

#Filter expression data for voxels that are in the pruned tree
treeMouseVoxels <- unlist(treeMouse$Get("voxels", filterFun = isLeaf))
dfExprMouse <- dfExprMouse[, (colnames(dfExprMouse) %in% treeMouseVoxels)]



# Build the human data tree ---------------------------------------------------

message("Importing human data tree...")

#Load human tree
#Samples are already mapped to nodes
load(str_c(pathHuman, "HumanExpressionTree_Samples.RData"))
treeHuman <- Clone(treeHumanExpr)
rm(treeHumanExpr)

#Remove white matter and ventricles
cutAtNodes <- c("white matter", "sulci & spaces")
pruneAnatTree(treeHuman, cutAtNodes, method = "AtNode")

#Filter expression data for samples that are in the pruned tree
treeHumanSamples <- treeHuman$`gray matter`$samples
dfExprHuman <- dfExprHuman[,colnames(dfExprHuman) %in% treeHumanSamples]



# Assign regional labels ------------------------------------------------------

message("Assigning labels...")

#Transpose data frames
dfExprMouse <- dfExprMouse %>% as.matrix() %>% t()
dfExprHuman <- dfExprHuman %>% as.matrix() %>% t()

#Assign gene names as columns. 
#Note: Necessary to do it this way if there are duplicate genes
colnames(dfExprMouse) <- genesMouse
colnames(dfExprHuman) <- genesHuman

#Convert back to data frame
#Note: This will create new names if there are duplicated genes. This is fine and desired.
dfExprMouse <- dfExprMouse %>% 
  as_tibble(rownames = "Voxels", .name_repair = "unique") %>% 
  column_to_rownames("Voxels")
dfExprHuman <- dfExprHuman %>% 
  as_tibble(rownames = "Samples", .name_repair = "unique") %>% 
  column_to_rownames("Samples")

#Load tree labels
load(str_c(pathHome, "Paper_Descriptive/Data/TreeLabels.RData"))

#Iterate over mouse parcellations
for (l in 1:length(listLabelsMouse)){
  
  #Prune tree to given level of aggregation
  treeMousePruned <- Clone(treeMouse)
  pruneAnatTree(treeMousePruned, listLabelsMouse[[l]], method = "BelowNode")
  
  #Get the region labels for every voxel
  dfExprMouse[,names(listLabelsMouse)[[l]]] <- labelRegions(rownames(dfExprMouse), treeMousePruned, "voxels")
}

#Iterate over human parcellations
for (l in 1:length(listLabelsHuman)){
  
  #Prune tree
  treeHumanPruned <- Clone(treeHuman)
  pruneAnatTree(treeHumanPruned, listLabelsHuman[[l]], method = "BelowNode")
  
  #Get the region labels
  dfExprHuman[,names(listLabelsHuman)[[l]]] <- labelRegions(rownames(dfExprHuman), treeHumanPruned, "samples")
}

#Remove identifiers from rownames
rownames(dfExprMouse) <- NULL
rownames(dfExprHuman) <- NULL


# Write to file ---------------------------------------------------------------

message("Writing to file...")

if (args[["savemouse"]] == "true") {
  outFileMouse <- str_c(str_remove(fileMouse, ".csv"), "_labelled", ".csv")  
  outPathMouse <- str_c(pathHome, "Paper_Descriptive/", "Data/", outFileMouse)  
  write_csv(dfExprMouse, path = outPathMouse)
}

if (args[["savehuman"]] == "true"){
  outFileHuman <- str_c(str_remove(fileHuman, ".csv"), "_labelled", ".csv")
  outPathHuman <- str_c(pathHome, "Paper_Descriptive/", "Data/", outFileHuman)
  write_csv(dfExprHuman, path = outPathHuman)
}
