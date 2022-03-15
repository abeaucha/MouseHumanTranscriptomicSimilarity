
library(optparse)

option_list <- list(
  make_option("--infile",
              type = "character")
)

args <- parse_args(OptionParser(option_list = option_list))

infile <- args[["infile"]]

Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
rmarkdown::render(infile, output_format = 'pdf_document')
