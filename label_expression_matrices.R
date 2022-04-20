# ----------------------------------------------------------------------------
# label_expression_matrices.R 
# Antoine Beauchamp
# 
# Label voxel-wise and sample-wise expression matrices
#
# Description
# -----------
# This script does two things:
# 1. It intersects the mouse and human gene sets with a set of homologous
#    genes
# 2. It labels the mouse and human voxel- and sample-wise expression matrices
#    with neuroanatomical labels from multiple sets of labels at different
#    levels of granularity in the neuroanatomical hierarchy. 


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(optparse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--mousematrix",
              type = "character",
              help = paste("Path to CSV file containing mouse voxel",
                           "expression matrix to label.")),
  make_option("--humanmatrix",
              type = "character",
              help = paste("Path to CSV file containing human sample",
                           "expression matrix to label.")),
  make_option("--mousetree",
              type = "character",
              help = paste("")),
  make_option("--humantree",
              type = "character",
              help = paste("")),
  make_option("--homologs",
              type = "character",
              help = paste("Path to CSV file containing mouse and human gene",
                           "homologs")),
  make_option("--outdir",
              type = "character",
              default = "data/",
              help = paste("[default %default]")),
  make_option("--savemouse",
              type = "character",
              default = "true",
              help = paste("Option to save labelled mouse data to file",
                           "[default %default]")),
  make_option("--savehuman",
              type = "character",
              default = "true",
              help = paste("Option to save labelled human data to file",
                           "[default %default]")),
  make_option("--verbose",
              default = "true",
              type = "character",
              help = "[default %default]")
)


#Parse command line args
args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args[["mousematrix"]])){
  stop("No input file given to --mousematrix")
}

if (is.null(args[["humanmatrix"]])){
  stop("No input file given to --humanmatrix")
}


if (!(args[["savemouse"]] %in% c("true", "false"))) {
  stop(paste("Argument --savemouse must be one of [true, false] (got ", 
             args[["savemouse"]], ")"))
}

if (!(args[["savehuman"]] %in% c("true", "false"))) {
  stop(paste("Argument --savehuman must be one of [true, false] (got ", 
             args[["savehuman"]], ")"))
}

if (!(args[['verbose']] %in% c('true', 'false'))) {
  stop(paste("Argument --verbose must be onoe of [true, false] (got ",
             args[['verbose']], ")"))
}


# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset("--file=") %>% 
  str_remove("--file=") %>% 
  dirname()

path_tree_tools <- file.path(working_dir, 
                             script_dir, 
                             "functions", 
                             "tree_tools.R")
source(path_tree_tools)

path_processing_tools <- file.path(working_dir,
                                   script_dir,
                                   "functions",
                                   "processing_tools.R")
source(path_processing_tools)


#' Label observations with neuroanatomical regions
#'
#' @param measurements (character vector) Unique identifiers for measurements
#' being labelled.
#' @param tree (data.tree tree) A tree whose leaf nodes are the labels you 
#' want to assign to the measurements.
#' @param treefield (character scalar) The name of the tree field that contains
#' the measurement identifiers for each node.
#'
#' @return (character vector) The neuroanatomical labels for each measurement.
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


# Main -----------------------------------------------------------------------

fileMouseMat <- args[["mousematrix"]]
fileMouseTree <- args[["mousetree"]]
fileHumanMat <- args[["humanmatrix"]]
fileHumanTree <- args[["humantree"]]
homologs <- args[["homologs"]]
verbose <- ifelse(args[["verbose"]] == 'true', TRUE, FALSE)

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

if (verbose) {message("Labelling data from mouse file: ", fileMouseMat)}
if (verbose) {message("Labelling data from human file: ", fileHumanMat)}

maskFlag <- str_extract(fileMouseMat, "mask[a-z]+")


# Importing and processing ---------------------------------------------------

if (verbose) {message("Importing data...")}

#Import expression matrices
dfExprMouse <- suppressMessages(data.table::fread(fileMouseMat, 
                                                  header = TRUE)) %>% 
  as_tibble()

dfExprHuman <- suppressMessages(data.table::fread(fileHumanMat,
                                                  header = TRUE)) %>% 
  as_tibble()

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


# Build the mouse data tree --------------------------------------------------

if (verbose) {message("Importing mouse data tree...")}

#Load DSURQE atlas and mask
dsurqe <- mincGetVolume("AMBA/data/imaging/DSURQE_CCFv3_labels_200um.mnc")

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

if (verbose) {message("Importing human data tree...")}

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

if (verbose) {message("Assigning labels...")}

#Transpose data frames
dfExprMouse <- dfExprMouse %>% as.matrix() %>% t()
dfExprHuman <- dfExprHuman %>% as.matrix() %>% t()

#Assign gene names as columns. 
#Note: Necessary to do it this way if there are duplicate genes
colnames(dfExprMouse) <- genesMouse
colnames(dfExprHuman) <- genesHuman

#Convert back to data frame
#Note: This will create new names if there are duplicated genes. 
#This is fine and desired.
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
  dfExprMouse[,names(listLabelsMouse)[[l]]] <- labelRegions(rownames(dfExprMouse),
                                                            treeMousePruned,
                                                            "voxels")
}

#Iterate over human parcellations
for (l in 1:length(listLabelsHuman)){
  
  #Prune tree
  treeHumanPruned <- Clone(treeHuman)
  pruneAnatTree(treeHumanPruned, listLabelsHuman[[l]], method = "BelowNode")
  
  #Get the region labels
  dfExprHuman[,names(listLabelsHuman)[[l]]] <- labelRegions(rownames(dfExprHuman),
                                                            treeHumanPruned, 
                                                            "samples")
}

#Remove identifiers from rownames
rownames(dfExprMouse) <- NULL
rownames(dfExprHuman) <- NULL


# Write to file --------------------------------------------------------------

if (verbose) {message("Writing to file...")}

if (args[["savemouse"]] == "true") {
  outFileMouse <- fileMouseMat %>% 
    basename() %>% 
    str_replace(".csv", "_labelled.csv")
  data.table::fwrite(dfExprMouse, file = str_c(args[['outdir']], outFileMouse))
}

if (args[["savehuman"]] == "true"){
  outFileHuman <- fileHumanMat %>% 
    basename() %>% 
    str_replace(".csv", "_labelled.csv")
  data.table::fwrite(dfExprHuman, file = str_c(args[['outdir']], outFileHuman))
}
