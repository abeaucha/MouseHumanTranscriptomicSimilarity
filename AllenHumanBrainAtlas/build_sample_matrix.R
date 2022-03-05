# BuildSampleExprMatrix.R ---------------------------------------------------------------
#
# This is a script to process the Allen Human Brain Atlas microarray gene expression 
# data and build a gene-by-sample expression matrix. The processing steps are loosely
# based on the paper "A practical guide to linking brain-wide gene expression and 
# neuroimaging data" by Aurina Arnatkeviciute, published in NeuroImage, 2019. 
# Two versions of the pipeline are available
#
# Output:
#
# Antoine Beauchamp
# Created: August 20th, 2021
# Edited: March 5th, 2022


# Libraries -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
library(optparse)


# Command line arguments ----------------------------------------------------------------

option_list <- list(
  make_option("--version",
              type = "numeric",
              default = 1,
              help = "Version of the data processing pipeline to run [default %default]"),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = "Verbose option for data processing pipeline [default %default]")
)

args <- parse_args(OptionParser(option_list = option_list))

if (!(args[["version"]] %in% c(1,2))){
  stop(str_c("Argument --version must be one of [1, 2] (got ", args[["version"]], ")"))
}

if (!(args[["verbose"]] %in% c("true","false"))){
  stop(str_c("Argument --verbose must be one of [true, false] (got ", args[["verbose"]], ")"))
}


# Paths ---------------------------------------------------------------------------------

#Set paths to AHBA directories
pathData <- "data/microarray/"
donorDirectories <- list.files(pathData)


# Functions -----------------------------------------------------------------------------

source("functions/processing_tools.R")


# Import --------------------------------------------------------------------------------

message("Importing data...")

#Import donor information and naming conventions
dfDonorInfo <- suppressMessages(read_csv("data/donors.csv"))) %>% 
  inner_join(tibble(donorFileID = str_remove(donorDirectories, "normalized_microarray_"),
                    donorPath = str_c(pathData, donorDirectories, "/")),
             by = "donorFileID") %>% 
  select(-donorFileID)
 
#Import data for each donor
listImport <- map2(dfDonorInfo$donorPath, dfDonorInfo$donorID, importDonorData)
names(listImport) <- dfDonorInfo$donorID


# Pipeline ------------------------------------------------------------------------------

#Apply processing pipeline to the data
listData <- processingPipeline(listImport,
                               version = args[["version"]],
                               verbose = ifelse(args[["verbose"]] == "true", TRUE, FALSE))


# Write ---------------------------------------------------------------------------------

message("Writing to file...")

#Write gene expression matrix
exprfile <- str_c("data/", "HumanExpressionMatrix_Samples_pipeline_v", args[["version"]], ".csv")
write_csv(x = listData[["GeneExpression"]],
          file = exprfile)

#Write sample information
samplefile <- str_c("data/", "SampleInformation_pipeline_v", args[["version"]], ".csv")
write_csv(x = listData[["SampleInfo"]],
          file = samplefile)
