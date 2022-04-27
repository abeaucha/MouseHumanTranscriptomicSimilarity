# ----------------------------------------------------------------------------
# extract_code.R
# Antoine Beauchamp
#
# Extract R code from .Rmd file

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--infile",
              type = "character")
)

args <- parse_args(OptionParser(option_list = option_list))

# Main -----------------------------------------------------------------------

infile <- args[["infile"]]

if (is.null(infile)){
  stop("No input file provided.") 
}

if (!stringr::str_detect(infile, ".Rmd$")) {
  stop("Input file must been an R Markdown file with extension .Rmd. Received: ", infile)
}

outfile <- stringr::str_replace(infile, ".Rmd", ".R")
knitr::purl(input = infile, output = outfile)