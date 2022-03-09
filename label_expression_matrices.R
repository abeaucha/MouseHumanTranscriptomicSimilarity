# label_data.R ------------------------------------------------------------
#
# 
# 
#
# Antoine Beauchamp
# Edited: March 9th, 2022

# Libraries -------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
library(RMINC)
library(data.tree)
library(optparse)


# Command line arguments ------------------------------------------------------

option_list <- list(
  make_option("--mousematrix",
              type = "character",
              help = "Path to CSV file containing mouse voxel expression matrix to label. [default %default]"),
  make_option("--humanmatrix",
              type = "character",
              help = "Path to CSV file containing human sample expression matrix to label. [default %default]"),
  make_option("--mousetree",
              type = "character"),
  make_option("--humantree",
              type = "character"),
  make_option("--homologs",
              type = "character",
              help = "Path to CSV file containing mouse and human gene homologs"),
  make_option("--outdir",
              type = "character",
              default = "data/"),
  make_option("--savemouse",
              type = "character",
              default = "true",
              help = "Option to save labelled mouse data to file [default %default]"),
  make_option("--savehuman",
              type = "character",
              default = "true",
              help = "Option to save labelled human data to file [default %default]")
)


# Functions -------------------------------------------------------------------

source("functions/tree_tools.R")
source("functions/processing_tools.R")

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

#Parse command line args
args = parse_args(OptionParser(option_list = option_list))

if (is.null(args[["mousematrix"]])){
  stop("No input file given to --mousematrix")
}

if (is.null(args[["humanmatrix"]])){
  stop("No input file given to --humanmatrix")
}


if (!(args[["savemouse"]] %in% c("true", "false"))) {
  stop(str_c("Argument --savemouse must be one of [true, false] (got ", args[["savemouse"]], ")"))
}

if (!(args[["savehuman"]] %in% c("true", "false"))) {
  stop(str_c("Argument --savehuman must be one of [true, false] (got ", args[["savehuman"]], ")"))
}


fileMouseMat <- args[["mousematrix"]]
fileMouseTree <- args[["mousetree"]]
fileHumanMat <- args[["humanmatrix"]]
fileHumanTree <- args[["humantree"]]
homologs <- args[["homologs"]]

if (!file.exists(fileMouseMat)) {
  stop("Mouse file ", fileMouseMat, " not found")
}

if (!file.exists(fileHumanMat)) {
  stop("Human file ", fileHumanMat, " not found")
}

if (!file.exists(fileMouseTree)) {
  stop("Mouse tree file ", fileMouseTree, " not found")
}

if (!file.exists(fileHumanTree)) {
  stop("Human tree file ", fileHumanTree, " not found")
}

if (!file.exists(homologs)) {
  stop("Homologs file ", homologs, " not found")
}

message("Labelling data from mouse file: ", fileMouseMat)
message("Labelling data from human file: ", fileHumanMat)

maskFlag <- str_extract(fileMouseMat, "mask[a-z]+")


# Importing and processing ----------------------------------------------------

message("Importing data...")

#Load expression data
dfExprMouse <- suppressMessages(read_csv(fileMouseMat))
dfExprHuman <- suppressMessages(read_csv(fileHumanMat))

#Subset genes for mouse-human homologs
listExpr <- intersectGeneHomologs(data = list(Mouse = dfExprMouse, 
                                              Human = dfExprHuman),
                                  homologs = homologs)

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
dsurqe <- mincGetVolume("AMBA/data/imaging/DSURQE_CCFv3_average_200um.mnc")

if (maskFlag == "masksagittal"){
  mask <- mincGetVolume("AMBA/data/imaging/sagittal_200um_coverage_bin0.8.mnc")
} else if (maskFlag == "maskcoronal") {
  mask <- mincGetVolume("AMBA/data/imaging/coronal_200um_coverage_bin0.8.mnc")
} else {
  stop(str_c("Invalid maskFlag value: ", maskFlag))
}

#Mask the DSURQE atlas
dsurqe_masked <- dsurqe[mask == 1]

#Load the mouse tree
load(fileMouseTree)

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
load(fileHumanTree)
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
load("data/TreeLabels.RData")

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
  outFileMouse <- fileMouseMat %>% 
    basename() %>% 
    str_replace(".csv", "_labelled.csv")
  write_csv(dfExprMouse, file = str_c(args[["outdir"]], outFileMouse))
}

if (args[["savehuman"]] == "true"){
  outFileHuman <- fileHumanMat %>% 
    basename() %>% 
    str_replace(".csv", "_labelled.csv")
  write_csv(dfExprHuman, file = str_c(args[["outdir"]], outFileHuman))
}
