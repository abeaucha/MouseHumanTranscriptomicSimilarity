# Directory: AllenHumanBrainAtlas

This directory contains scripts and data related to the Allen Human Brain Atlas from the Allen Institute for Brain Science. 

## Sub-directories

**AtlasMappings**: Directory containing mappings from the AHBA microarray samples to various atlases.

**Data**: Contains data related to human gene expression.

**Functions**: Functions

## Files

**BuildHumanExprTree.R**: R script that builds a hierarchical anatomical tree containing gene expression data at every node. This can be used to build a tree for individual donors, or for the donor-combined data. The choice of atlas can also be specified. Trees are saved in the Data sub-directory.

**BuildSampleExprMatrix.py**: Python script that builds a gene-by-sample expression matrix using the data from all of the donors. This matrix is saved in the Data sub-directory.
ls


