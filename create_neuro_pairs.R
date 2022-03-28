# ----------------------------------------------------------------------------
# create_neuro_pairs.R
# Author: Antoine Beauchamp
# Created: August 4th, 2021
#
# Description
# -----------
# This script defines a set of canonical neuroanatomical pairs between the
# mouse and human brain atlases

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))

# Command line arguments -----------------------------------------------------

option_list(
  make_option("--outdir",
              type = "character",
              default = "data/",
              help = paste("Directory in which to save the CSV file.",
                           "[default %default]")),
  make_option("--outfile",
              type = "character",
              default = "MouseHumanMatches_H88M67.csv",
              help = paste("Name of CSV file containing mouse-human",
                           "neuroanatomical pairs. [default %default]"))
)

args <- parse_args(OptionParser(option_list = option_list))

# Mouse labels 67 / Human labels 88 -------------------------------------------

dfMouseHumanMatches_H88M67 <- tibble::tibble(Mouse = c("Claustrum", 
                                                       "Piriform area",
                                                       "Subiculum",
                                                       "Field CA1", 
                                                       "Field CA2",
                                                       "Field CA3",
                                                       "Dentate gyrus",
                                                       "Anterior cingulate area",
                                                       "Primary auditory area",
                                                       "Primary motor area",
                                                       "Primary somatosensory area",
                                                       "Visual areas",
                                                       "Pallidum",
                                                       "Striatum ventral region",
                                                       "Caudoputamen",
                                                       "Cortical subplate-other",
                                                       "Inferior colliculus",
                                                       "Superior colliculus, sensory related",
                                                       "Medulla",
                                                       "Pons",
                                                       "Hypothalamus",
                                                       "Thalamus",
                                                       "Lingula (I)",
                                                       "Declive (VI)",
                                                       "Folium-tuber vermis (VII)",
                                                       "Pyramus (VIII)",
                                                       "Uvula (IX)", 
                                                       "Nodulus (X)", 
                                                       "Simple lobule",
                                                       "Crus 1",
                                                       "Crus 2",
                                                       "Paramedian lobule",
                                                       "Copula pyramidis",
                                                       "Flocculus",
                                                       "Paraflocculus",
                                                       "Cerebellar nuclei"),
                                             Human = c("claustrum",
                                                       "piriform cortex",
                                                       "subiculum",
                                                       "CA1 field",
                                                       "CA2 field",
                                                       "CA3 field",
                                                       "dentate gyrus",
                                                       "cingulate gyrus",
                                                       "Heschl's gyrus", 
                                                       "precentral gyrus", 
                                                       "postcentral gyrus",
                                                       "cuneus",
                                                       "globus pallidus",
                                                       "nucleus accumbens",
                                                       "caudate nucleus",
                                                       "amygdala",
                                                       "inferior colliculus",
                                                       "superior colliculus",
                                                       "myelencephalon",
                                                       "pons",
                                                       "hypothalamus",
                                                       "thalamus",
                                                       "vermal I-II",
                                                       "vermal VI", 
                                                       "vermal VIIAf",
                                                       "vermal VIIIA", 
                                                       "vermal IX", 
                                                       "vermal X",
                                                       "VI",
                                                       "crus I",
                                                       "crus II",
                                                       "VIIB",
                                                       "VIIIA", 
                                                       "X",
                                                       "IX",
                                                       "cerebellar nuclei"))


outfile <- file.path(args[["outdir"]], args[["outfile"]])
readr::write_csv(x = dfMouseHumanMatches_H88M67, 
                 file = outfile)

