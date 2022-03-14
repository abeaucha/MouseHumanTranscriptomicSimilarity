# process_labelled_matrices.R
#
# Antoine Beauchamp
# Created: August 25th, 2021
# Edited: March 9th, 2022

#Libraries
suppressPackageStartupMessages(library(tidyverse))
library(optparse)

# Command line arguments 

option_list <- list(
  make_option("--infile",
              type = "character"),
  make_option("--scale",
              type = "character",
              default = "true"),
  make_option("--aggregate",
              type = "character",
              default = "false"),
  make_option("--nlabels",
              type = "integer"),
  make_option("--outdir",
              default = "data/",
              type = "character")
)


#Functions
source("functions/processing_tools.R")

args <- parse_args(OptionParser(option_list = option_list))

if (!(args[["scale"]] %in% c("true", "false"))) {
  stop()
}

if (!(args[["aggregate"]] %in% c("true", "false"))) {
  stop()
}

if (args[["scale"]] == "false" & args[["aggregate"]] == "false"){
  stop()
}


#Import data
dfExpression <- suppressMessages(read_csv(args["infile"]))

#Extract genes list from data
genes <- colnames(dfExpression)[!str_detect(colnames(dfExpression), "Region")]

if (args[["scale"]] == "true") {
  
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

write_csv(dfExpression,
          file = str_c(args[["outdir"]],outfile))
