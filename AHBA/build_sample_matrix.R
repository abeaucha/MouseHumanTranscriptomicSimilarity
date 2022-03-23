# ----------------------------------------------------------------------------
# build_sample_matrix.R
# Author: Antoine Beauchamp
# Created: August 20th, 2021
#
# Description
# -----------
# This is a script to process the Allen Human Brain Atlas microarray gene 
# expression data and build a gene-by-sample expression matrix. 
# The processing steps are loosely based on the paper 
# "A practical guide to linking brain-wide gene expression and 
# neuroimaging data" by Aurina Arnatkeviciute, published in NeuroImage, 2019. 
# Two versions of the pipeline are available.


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--version",
              type = "numeric",
              default = 1,
              help = paste("Version of the data processing pipeline to run",
                           "[default %default]")),
  make_option("--verbose",
              type = "character",
              default = "true",
              help = paste("Verbose option for data processing pipeline",
                           "[default %default]"))
)

args <- parse_args(OptionParser(option_list = option_list))

if (!(args[["version"]] %in% c(1,2))){
  stop(paste("Argument --version must be one of [1, 2]",
             "(got", args[["version"]], ")"))
}

if (!(args[["verbose"]] %in% c("true","false"))){
  stop(paste("Argument --verbose must be one of [true, false]",
             "(got", args[["verbose"]], ")"))
}


# Functions ------------------------------------------------------------------

source("processing_tools_AHBA.R")


# Import ---------------------------------------------------------------------

#Set paths to AHBA directories
pathData <- "data/microarray/"
donorDirectories <- list.files(pathData)

message("Importing data...")

#Import donor information and naming conventions
dfDonorInfo <- suppressMessages(read_csv("data/donors.csv")) %>% 
  inner_join(tibble(donorFileID = str_remove(donorDirectories, 
                                             "normalized_microarray_"),
                    donorPath = str_c(pathData, 
                                      donorDirectories, 
                                      "/")),
             by = "donorFileID") %>% 
  select(-donorFileID)

#Import data for each donor
listImport <- map2(dfDonorInfo$donorPath, 
                   dfDonorInfo$donorID,
                   importDonorData)
names(listImport) <- dfDonorInfo$donorID


# Pipeline -------------------------------------------------------------------

#Apply processing pipeline to the data
listData <- processingPipeline(listImport,
                               version = args[["version"]],
                               verbose = ifelse(args[["verbose"]] == "true", 
                                                TRUE,
                                                FALSE))


# Write ----------------------------------------------------------------------

message("Writing to file...")

#Write gene expression matrix
exprfile <- paste0("HumanExpressionMatrix_samples_pipeline_v", args[["version"]], ".csv")
write_csv(x = listData[["GeneExpression"]],
          file = paste0("data/", exprfile))

#Write sample information
samplefile <- paste0("SampleInformation_pipeline_v", args[["version"]], ".csv")
write_csv(x = listData[["SampleInfo"]],
          file = paste0("data/", samplefile))
