#!/bin/bash

# ----------------------------------------------------------------------------
# pipeline_main
# Author: Antoine Beauchamp
#
# Pipeline to download and process all data necessary to generate the paper.

#On MICe machine
module purge

source activate_venv.sh

#Download and process AMBA data
source build_AMBA_data.sh

#Download and process AHBA data
source build_AHBA_data.sh

#Jointly process AMBA and AHBA expression matrices
source process_expression_data.sh

#Generate 500 gene expression latent spaces
source generate_latent_spaces.sh

deactivate
