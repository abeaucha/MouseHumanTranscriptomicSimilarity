#!/bin/bash

# ----------------------------------------------------------------------------
# build_AMBA_data
# Author: Antoine Beauchamp
#
# Pipeline to generate necessary data files from the Allen Mouse Brain Atlas
#
# Pipeline steps:
# 1. Resample DSURQE imaging files from its common space to CCFv3. 
# 2. Download the AMBA coronal in-situ hybridization data set from the web
# 3. Download the AMBA sagittal in-situ hybridization data set from the web
# 4. Build a gene-by-voxel expression matrix using the coronal data set with 
#    a bilateral coronal imaging mask
# 5. Build a gene-by-voxel expression matrix using the coronal data set with 
#    a unilateral sagittal imaging mask
# 6. Build a gene-by-voxel expression matrix using the sagittal data set with 
#    a unilateral sagittal imaging mask
# 7. Build a gene-by-region expression matrix from the coronal voxel-wise 
#    expression matrix using the DSURQE atlas.
# 8. Build a gene expression tree by combining the coronal voxel-wise
#    expression matrix with the AMBA hierarchical ontology. 

# On MICe machines
#module purge
#source AMBA/resample_DSURQE_CCFv3

# Activate the python virtual environment
source activate_venv

# Download AHBA data from the web
echo "Downloading AMBA coronal in-situ hybridization data set..."
python3 AMBA/download_AMBA.py \
	--dataset coronal \
	--outdir AMBA/data/expression/ \
	--metadata AMBA_metadata_coronal.csv \
	--parallel true \
	--nproc 12

echo "Downloading AMBA sagittal in-situ hybridization data set..."
python3 AMBA/download_AMBA.py \
	--dataset sagittal \
	--outdir AMBA/data/expression/ \
	--metadata AMBA_metadata_sagittal.csv \
	--parallel true \
	--nproc 12

# Build voxel expression matrix 
echo "Building mouse gene-by-voxel expression matrix using coronal data with coronal mask..."
python3 AMBA/build_voxel_matrix.py \
	--datadir AMBA/data/expression/ \
	--outdir AMBA/data/ \
	--imgdir AMBA/data/imaging/ \
	--dataset coronal \
	--mask coronal \
	--log2 true \
	--groupexp true \
	--threshold 0.2 \
	--impute true \
	--parallel true \
	--nproc 4 

echo "Building mouse gene-by-voxel expression matrix using coronal data with sagittal mask..."
python3 AMBA/build_voxel_matrix.py \
	--datadir AMBA/data/expression/ \
	--outdir AMBA/data/ \
	--imgdir AMBA/data/imaging/ \
	--dataset coronal \
	--mask sagittal \
	--log2 true \
	--groupexp false \
	--threshold 0.2 \
	--impute true \
	--parallel true \
	--nproc 4 

echo "Building mouse gene-by-voxel expression matrix using sagittal data with sagittal mask..."
python3 AMBA/build_voxel_matrix.py \
	--datadir AMBA/data/expression/ \
	--outdir AMBA/data/ \
	--imgdir AMBA/data/imaging/ \
	--dataset sagittal \
	--mask sagittal \
	--log2 true \
	--groupexp true \
	--threshold 0.2 \
	--impute true \
	--parallel true \
	--nproc 4 

# Build regional expression matrix using DSURQE
echo "Building mouse gene-by-region expression matrix using the DSURQE atlas..."
python3 AMBA/build_region_matrix.py \
	--datadir AMBA/data/ \
	--imgdir AMBA/data/imaging \
	--infile MouseExpressionMatrix_voxel_coronal_maskcoronal_log2_grouped_imputed.csv \
	--outfile MouseExpressionMatrix_ROI_DSURQE_coronal_maskcoronal_log2_grouped_imputed.csv \
	--mask coronal_200um_coverage_bin0.8.mnc \
	--labels DSURQE_CCFv3_labels_200um.mnc \
	--defs DSURQE_40micron_R_mapping_long.csv
 
# Build mouse data tree
echo "Building mouse neuroanatomical tree..."
Rscript AMBA/build_mouse_tree.R \
	--datadir AMBA/data/ \
	--imgdir AMBA/data/imaging/ \
	--infile MouseExpressionMatrix_ROI_DSURQE_coronal_maskcoronal_log2_grouped_imputed.csv \
	--outfile MouseExpressionTree_DSURQE.RData \
	--treefile DSURQE_tree.json \
	--defs DSURQE_40micron_R_mapping_long.csv

deactivate
