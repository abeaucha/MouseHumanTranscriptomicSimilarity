# ----------------------------------------------------------------------------
# render_pdf.R
# Antoine Beauchamp
#
# Render .pdf file from .Rmd file

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))

# Command line arguments -----------------------------------------------------

option_list <- list(
  make_option("--infile",
              type = "character")
)

args <- parse_args(OptionParser(option_list = option_list))

# Main -----------------------------------------------------------------------

Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
rmarkdown::render(args[["infile"]], output_format = 'pdf_document')
