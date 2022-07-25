#!/bin/bash

# ----------------------------------------------------------------------------
# process_cortical_expression_data
# Author: Antoine Beauchamp
#
# Process mouse and human cortical expression data sets.

#Data directory
datadir=data/isocortex/

#Create a mouse voxel-wise expression matrix with only isocortical labels
Rscript subset_cortical_labels.R

# Normalize and aggregate cortical data

echo "Normalizing mouse isocortical expression matrix..."
Rscript process_labelled_matrix.R \
	--infile ${datadir}MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled.csv \
	--scale true \
	--aggregate false \
	--outdir ${datadir} \
	--verbose true

echo "Normalizing human isocortical expression matrix..."
Rscript process_labelled_matrix.R \
	--infile ${datadir}HumanExpressionMatrix_samples_pipeline_abagen_labelled.csv \
	--scale true \
	--aggregate false \
	--outdir ${datadir} \
	--verbose true

echo "Normalizing and aggregating mouse isocortical expression matrix..."
Rscript process_labelled_matrix.R \
	--infile ${datadir}MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed_labelled.csv \
	--scale true \
	--aggregate true \
	--nlabels 67 \
	--outdir ${datadir} \
	--verbose true

echo "Normalizing and aggregating human isocortical expression matrix..."
Rscript process_labelled_matrix.R \
	--infile ${datadir}HumanExpressionMatrix_samples_pipeline_abagen_labelled.csv \
	--scale true \
	--aggregate true \
	--nlabels 88 \
	--outdir ${datadir} \
	--verbose true