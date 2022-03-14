

library(optparse)
library(stringr)

option_list <- list(
  make_option("--infile",
              type = "character")
)

args <- parse_args(OptionParser(option_list = option_list))

infile <- args[["infile"]]

if (is.null(infile)){
  stop("No input file provided") 
}

if (!str_detect(infile, ".Rmd$")) {
  stop("Input file must been an R Markdown file with extension .Rmd. Received: ", infile)
}

outfile <- str_replace(infile, ".Rmd", ".R")
knitr::purl(input = infile, output = outfile)