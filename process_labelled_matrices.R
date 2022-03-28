# ----------------------------------------------------------------------------
# process_labelled_matrices.R
# Antoine Beauchamp
# Created: August 25th, 2021
#
# Description
# -----------

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--infile",
              type = "character",
              help = paste("Path to CSV file containing labelled expression",
                           "matrix.")),
  make_option("--scale",
              type = "character",
              default = "true",
              help = paste("Option to scale expression data",
                           "[default %default]")),
  make_option("--aggregate",
              type = "character",
              default = "false",
              help = paste("Option to aggregate expression data under a",
                           "set of atlas labels [default %default]")),
  make_option("--nlabels",
              type = "integer",
              help = paste("Number of labels in the atlas used to",
                           "aggregate the data. Ignored if --aggregate",
                           "is false.")),
  make_option("--outdir",
              default = "data/",
              type = "character",
              help = paste("Directory in which to save processed data.",
                           "[default %default]")),
  make_option("--verbose",
              default = "true",
              type = "character",
              help = "[default %default]")
)

args <- parse_args(OptionParser(option_list = option_list))

if (!(args[["scale"]] %in% c("true", "false"))) {
  stop()
}

if (!(args[["aggregate"]] %in% c("true", "false"))) {
  stop()
}

if (args[["scale"]] == "false" & args[["aggregate"]] == "false"){
  stop(paste("Both --scale and --aggregate were false. This does nothing."))
}

# Functions ------------------------------------------------------------------

working_dir <- getwd()

script_dir <- commandArgs() %>% 
  str_subset("--file=") %>% 
  str_remove("--file=") %>% 
  dirname()

path_processing_tools <- file.path(working_dir,
                                   script_dir,
                                   "functions",
                                   "processing_tools.R")
source(path_processing_tools)


# Main -----------------------------------------------------------------------

verbose <- ifelse(args[["verbose"]] == 'true', TRUE, FALSE)

if(verbose){message("Importing data...")}

#Import data
dfExpression <- suppressMessages(data.table::fread(args[["infile"]],
                                                   header = TRUE)) %>% 
  as_tibble()

#Extract genes list from data
genes <- colnames(dfExpression)[!str_detect(colnames(dfExpression), "Region")]

if (args[["scale"]] == "true") {
  
  if(verbose){message("Scaling data...")}
  
  #Extract labels from data frame
  dfLabels <- dfExpression %>% select(contains("Region"))
  
  #Normalize data
  dfExpression <- dfExpression %>% 
    select(-contains("Region")) %>% 
    as.matrix() %>% 
    scaler(axis = "rows") %>% 
    scaler(scale = FALSE, axis = "columns") %>% 
    as_tibble() %>% 
    bind_cols(dfLabels)
  
  outfile <- args[["infile"]] %>% 
    basename() %>% 
    str_replace(".csv", "_scaled.csv")
  
}

if (args[["aggregate"]] == "true") {
  
  if(verbose){message("Aggregating data...")}
  
  labels <- str_c("Region", args[["nlabels"]])
  
  if (!(labels %in% colnames(dfExpression))) {
    stop()
  }
  
  #Aggregate mouse expression data under label set
  dfExpression <- dfExpression %>% 
    select(Region = all_of(labels), all_of(genes)) %>% 
    group_by(Region) %>% 
    summarise_all(mean) %>% 
    ungroup()
  
  outfile <- args[["infile"]] %>% 
    basename() %>% 
    str_extract("^[a-zA-Z]*") %>% 
    str_c("_ROI_", labels, ".csv")
  
  if (args[["scale"]] == "true") {
    outfile <- str_replace(outfile, ".csv", "_scaled.csv")
  }
  
}

if(verbose){message("Writing to file...")}

write_csv(dfExpression,
          file = file.path(args[["outdir"]], outfile))
