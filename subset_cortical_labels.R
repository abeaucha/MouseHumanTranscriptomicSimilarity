#Packages
suppressPackageStartupMessages(library(tidyverse))

#Output directory
outdir='data/isocortex/'

#Create outdir if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

#Input and output mouse files
infile <- 'data/MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled.csv'
outfile <- file.path(outdir, basename(infile))

#Filter mouse input matrix for isocortical labels and write out
data.table::fread(infile, header = TRUE) %>% 
  as_tibble() %>% 
  filter(Region11 == 'Isocortex') %>% 
  data.table::fwrite(outfile)

#Input and output human files
infile <- 'data/HumanExpressionMatrix_samples_pipeline_abagen_labelled.csv'
outfile <- file.path(outdir, basename(infile))

#Filter human input matrix for cortical labels and write out
data.table::fread(infile, header = TRUE) %>% 
  as_tibble() %>% 
  filter(Region16 %in% c("limbic lobe", 
                         "frontal lobe", 
                         "insula", 
                         "occipital lobe",
                         "parietal lobe", 
                         "temporal lobe")) %>% 
  filter(!(Region88 %in% c("claustrum",
                           "dentate gyrus",
                           "CA1 field",
                           "CA2 field",
                           "CA3 field",
                           "CA4 field",
                           "subiculum"))) %>% 
  data.table::fwrite(outfile)

