# build_region_matrix.py -------------------------------------
#
#
#
# Antoine Beauchamp
# Created: January 31st, 2022
# Edited: March 7th, 2022
# --------------------------------------------------------------

# Packages -----------------------------------------------------

import argparse
import numpy as np
import pandas as pd
from pyminc.volumes.factory import *


# Functions ---------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--dataset',
        type = str,
        default = 'coronal',
        choices = ['coronal', 'sagittal'],
        help = 'AMBA dataset to import'
    )
    
    args = vars(parser.parse_args())
    
    return args


def main():
    
    #Load command line arguments
    args = parse_args()
    
    #Flags
    dataset = args['dataset']
    
    print("Using {} dataset".format(dataset))
          
    
    # Import expression data --------------------------------------
    
    #File containing voxel expression matrix
    infile = 'MouseExpressionMatrix_voxel_'+dataset+'_mask'+dataset+'_imputed.csv'
    
    print("Importing voxel expression matrix: {} ...".format(infile))
    
    #Import voxel expression matrix
    dfExprVoxel = pd.read_csv('data/{}'.format(infile))
    
    #Extract gene names from data frame
    genes = dfExprVoxel['Gene']

    #Extract voxels from data frame and convert to numpy array
    npExprVoxel = np.array(dfExprVoxel.loc[:,dfExprVoxel.columns != 'Gene'])
    
    #Transpose numpy array so that voxels are rows
    npExprVoxel = np.transpose(npExprVoxel)

    #Convert transposed numpy array to data frame. Set genes to be column names.
    dfExprVoxelTranspose = pd.DataFrame(npExprVoxel, columns=genes)
          
        
        
    # Import imaging data --------------------------------------
    
    #Import image mask, flatten and convert to numpy array
    if dataset == 'sagittal':
        maskfile = 'sagittal_200um_coverage_bin0.8.mnc'
    else:
        maskfile = 'coronal_200um_coverage_bin0.8.mnc'
    
    print("Importing mask volume: {} ...".format(maskfile))
    
    maskVol = volumeFromFile('data/imaging/{}'.format(maskfile))
    maskArray = np.array(maskVol.data.flatten())
    maskVol.closeVolume()
    
    
    #File containing atlas labels
    labelfile = 'DSURQE_CCFv3_labels_200um.mnc'
    
    print("Importing DSURQE label volume: {} ...".format(labelfile))
    
    #Import DSURQE label volume, flatten, and convert to numpy array
    labelVol = volumeFromFile('data/imaging/{}'.format(labelfile))
    labelArray = np.array(labelVol.data.flatten())
    labelVol.closeVolume()

    #Mask the label array
    labelArrayMasked = labelArray[maskArray == 1]

    
    #File containing atlas definitions
    defsfile = 'DSURQE_40micron_R_mapping_long.csv'
    
    print("Importing DSURQE label definitions: {} ...".format(defsfile))
    
    #Import DSURQE label definitions
    dfAtlasDefs = pd.read_csv('data/imaging/{}'.format(defsfile))
    
    
    
    # Match labels to voxels -----------------------------------------
        
    print("Matching atlas labels to voxels...")

    #Initialize empty string array to match DSURQE regions to voxels
    structArray = np.empty(len(labelArrayMasked), dtype = 'U100')

    #Iterate over DSURQE regions and match to voxels
    for i, row in dfAtlasDefs.iterrows():
        lab = row['Label']
        indlab = labelArrayMasked == row['Label']
        
        if lab == 0:
            structArray[indlab] = np.nan
        else:
            structArray[indlab] = row['Structure']
    
    #Assign DSURQE regions to new column in the data frame
    dfExprVoxelTranspose['Region'] = structArray
    
    #Remove voxels that are empty
    dfExprVoxelTranspose = dfExprVoxelTranspose.loc[dfExprVoxelTranspose['Region'] != 'nan']
    
    
    
    # Aggregate expression data -------------------------------------------
    
    print("Aggregating expression data...")
    
    #Aggregate expression data by ROI
    dfExprRegionTranspose = dfExprVoxelTranspose.groupby('Region').aggregate(np.mean)

    #Transpose data frame so that rows are genes
    dfExprRegion = np.transpose(dfExprRegionTranspose)

    #Gene names are stored in the index
    #Label the index as Gene
    dfExprRegion.index.name = 'Gene'
    

    
    # Write -------------------------------------------
    
    print("Writing to file...")
    
    outfile = 'MouseExpressionMatrix_ROI_DSURQE_'+dataset+'_mask'+dataset+'.csv'
    
    #Write regional expression matrix to csv
    dfExprRegion.to_csv('data/{}'.format(outfile))
    
    return

    
    
if __name__ == '__main__':
    main()