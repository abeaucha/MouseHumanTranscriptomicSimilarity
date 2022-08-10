#!/bin/bash

# ----------------------------------------------------------------------------
# pipeline_main
# Author: Antoine Beauchamp
#
# Pipeline to download and process all data necessary for the analysis and
# manuscript.

#On MICe machine
module purge

source activate_venv

#Download and process AMBA data
source build_AMBA_data

#Download and process AHBA data
source build_AHBA_data

#Jointly process AMBA and AHBA expression matrices
source process_expression_data

#Generate 500 gene expression latent spaces
source generate_latent_spaces

deactivate
