# Directory: AllenGeneExpression_Mouse

This directory contains data related to the Allen Mouse Brain Atlas (AMBA).

## Sub-directories

**MRI**: Directory containing MRI data associated with the AMBA

## Files

**ABIMouseGeneExp_Coronal.txt**: List of MINC file names in the coronal ISH experiment.

**ABIMouseGeneExp_Sagittal.txt**: List of MINC file names in the sagittal ISH experiment.

**MouseExpressionMatrix_DSURQE.csv**: Gene-by-region expression matrix with regions drawn from the DSURQE atlas. Expressions haven't been normalized in any way.

**MouseExpressionMatrix_Voxel.csv**: Gene-by-voxel expression matrix. Un-normalized. This matrix was built using the coronal_200um_coverage_bin0.8.mnc mask in the MRI sub-directory. 

**MouseExpressionTree_DSURQE.RData**: Hierarchical anatomy tree containing gene expression values at every node.

**VoxelExpressionMatrix_YohanYee.npz**: (Symlink). Voxel expression matrix built by Yohan Yee. 

**VoxelExpressionmatrix_Metadata_YohanYee.npz**: (Symlink). Metadata associated with Yohan's voxel expression matrix.


