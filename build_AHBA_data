#!/bin/bash

# ----------------------------------------------------------------------------
# build_AHBA_data
# Author: Antoine Beauchamp
#
# Pipeline to generate necessary data files from the Allen Human Brain Atlas
#
# Pipeline steps:
# 1. Download the AHBA ontology and microarray data sets from the web
# 2. Build the gene-by-sample expression matrix using data from all donors
# 3. Build a gene expression tree using the sample expression matrix

# On MICe machines
module purge

# Activate the python virtual environment
source activate_venv

# Download AHBA data from the web
echo "Downloading AHBA data..."
python3 AHBA/download_AHBA.py \
	--outdir AHBA/data/

# Build sample expression matrix
echo "Building human gene-by-sample expression matrix..." 
python3 AHBA/build_sample_matrix.py \
	--datadir AHBA/data/microarray/ \
	--outdir AHBA/data/ \
	--donorsfile AHBA/data/donors.csv \
	--ibf-threshold 0.5 \
	--probe-selection diff_stability \
	--donor-probes aggregate \
	--sample-norm srs \
	--gene-norm srs \
	--verbose true
    
# Build expression tree
echo "Building human neuroanatomical tree..."
Rscript AHBA/build_human_tree.R \
	--datadir AHBA/data/ \
	--infile HumanExpressionMatrix_samples_pipeline_abagen.csv \
	--samplefile SampleInformation_pipeline_abagen.csv \
	--treefile AHBA_hierarchy_definitions.json \
	--outfile HumanExpressionTree.RData

deactivate
