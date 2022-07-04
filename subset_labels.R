library(tidyverse)

infile <- 'data/MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled_scaled.csv'
outfile <- 'data/isocortex/MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled_scaled.csv'

df_mouse <- data.table::fread(infile, header = TRUE) %>% 
  as_tibble() %>% 
  filter(Region11 == 'Isocortex') %>% 
  data.table::fwrite(outfile)

