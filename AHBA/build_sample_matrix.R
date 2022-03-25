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

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))


# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--microarraydir",
              type = "character",
              default = "data/microarray/",
              help = paste("Directory containing the microarray data sets",
                           "from the Allen Human Brain Atlas", 
                           "[default %default]")),
  make_option("--datadir",
              type = "character",
              default = "data/", 
              help = paste("Directory containing additional data files",
                           "This also acts as the output directory",
                           "[default %default]")),
  make_option("--donorsfile",
              type = "character",
              default = "donors.csv",
              help = paste("CSV file containing donor name conventions",
                           "[default %default]")),
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

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset("--file=") %>% 
  str_remove("--file=") %>% 
  dirname()

path_tools <- str_c(working_dir, script_dir, "processing_tools_AHBA.R",
                    sep = "/")

source(path_tools)


# Import ---------------------------------------------------------------------

dirMicroarray <- args[["microarraydir"]]
dirData <- args[["datadir"]]
fileDonors <- args[["donorsfile"]]

#Set paths to AHBA directories
donorDirectories <- list.files(dirMicroarray)

message("Importing data...")

#Import donor information and naming conventions
dfDonorInfo <- suppressMessages(read_csv(str_c(dirData, fileDonors))) %>% 
  inner_join(tibble(donorFileID = str_remove(donorDirectories, 
                                             "normalized_microarray_"),
                    donorPath = str_c(dirMicroarray, 
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
          file = paste0(dirData, exprfile))

#Write sample information
samplefile <- paste0("SampleInformation_pipeline_v", args[["version"]], ".csv")
write_csv(x = listData[["SampleInfo"]],
          file = paste0(dirData, samplefile))
