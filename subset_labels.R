suppressPackageStartupMessages(library(tidyverse))

infile <- 'data/MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled_scaled.csv'
outfile <- 'data/isocortex/MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled_scaled.csv'

data.table::fread(infile, header = TRUE) %>% 
  as_tibble() %>% 
  filter(Region11 == 'Isocortex') %>% 
  data.table::fwrite(outfile)

infile <- 'data/MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled.csv'
outfile <- 'data/isocortex_scaled/MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled.csv'

data.table::fread(infile, header = TRUE) %>% 
  as_tibble() %>% 
  filter(Region11 == 'Isocortex') %>% 
  data.table::fwrite(outfile)

infile <- 'data/HumanExpressionMatrix_samples_pipeline_abagen_labelled.csv'
outfile <- 'data/isocortex_scaled/HumanExpressionMatrix_samples_pipeline_abagen_labelled.csv'

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

