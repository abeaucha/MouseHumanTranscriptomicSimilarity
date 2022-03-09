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
import os
import numpy as np
import pandas as pd
from pyminc.volumes.factory import *


# Functions ---------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/',
        help = 'Directory containing expression matrices'
    )
    
    parser.add_argument(
        '--infile',
        type = str,
        help = 'CSV file containing expression matrix'
    )
    
    parser.add_argument(
        '--mask',
        type = str
    )
    
    parser.add_argument(
        '--labels',
        type = str
    )
    
    parser.add_argument(
        '--defs',
        type = str
    )
    
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
    datadir = args['datadir']
    infile = args['infile']
    mask = args['mask']
    labels = args['labels']
    defs = args['defs']
#     dataset = args['dataset']
    
#     print("Using {} dataset".format(dataset))
          
    
    # Import expression data --------------------------------------
    
    #File containing voxel expression matrix
#     infile = 'MouseExpressionMatrix_voxel_'+dataset+'_mask'+dataset+'_imputed.csv'
    
    print("Importing gene-by-voxel expression matrix: {} ...".format(infile))
    
    #Import voxel expression matrix
    dfExprVoxel = pd.read_csv(os.path.join(datadir, infile))
    
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
    print("Importing mask volume: {} ...".format(mask))

    maskVol = volumeFromFile(os.path.join('data', 'imaging', mask))
    maskArray = np.array(maskVol.data.flatten())
    maskVol.closeVolume()
    
    
    #File containing atlas labels
#     labelfile = 'DSURQE_CCFv3_labels_200um.mnc'
    
    print("Importing DSURQE label volume: {} ...".format(labels))
    
    #Import DSURQE label volume, flatten, and convert to numpy array
    labelVol = volumeFromFile(os.path.join('data', 'imaging', labels))
    labelArray = np.array(labelVol.data.flatten())
    labelVol.closeVolume()

    #Mask the label array
    labelArrayMasked = labelArray[maskArray == 1]

    
    #File containing atlas definitions
#     defsfile = 'DSURQE_40micron_R_mapping_long.csv'
    
    print("Importing DSURQE label definitions: {} ...".format(defs))
    
    #Import DSURQE label definitions
    dfAtlasDefs = pd.read_csv(os.path.join('data', 'imaging', defs))
    
    
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
    dfExprRegion.to_csv(os.path.join(datadir, outfile))
    
    return

    
    
if __name__ == '__main__':
    main()