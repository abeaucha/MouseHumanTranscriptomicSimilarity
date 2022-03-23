# ----------------------------------------------------------------------------
# build_region_matrix.py
# Author: Antoine Beauchamp
# Created: January 31st, 2022

"""
Build a gene-by-region expression matrix

Description
-----------
This script imports a gene-by-voxel expression matrix from a CSV file
and aggregates the expression values for every region in an atlas. 
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import numpy as np
import pandas as pd
from pyminc.volumes.factory import *

# Command line arguments -----------------------------------------------------

def parse_args():
   
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/',
        help = "Directory containing expression matrix CSV files"
    )
    
    parser.add_argument(
        '--infile',
        type = str,
        help = "Name of CSV file containing the voxel-wise expression matrix"
    )
    
    parser.add_argument(
        '--outfile',
        type = str,
        help = "Name of CSV file in which to save regional expression matrix"
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        help = ("Name of MINC file containing the mask used to build the "
                "voxel-wise expression matrix provided to --infile")
    )
    
    parser.add_argument(
        '--labels',
        type = str,
        help = ("Name of MINC file containing the atlas labels to use. "
                "Must be in CCFv3 space.")
    )
    
    parser.add_argument(
        '--defs',
        type = str,
        help = ("Name of CSV file containing the names of the "
                "neuroanatomical regions corresponding to the "
                "atlas labels in --labels")
    )
    
    args = vars(parser.parse_args())
    
    return args

# Main -----------------------------------------------------------------------

def main():
    
    #Load command line arguments
    args = parse_args()
    
    #Flags
    datadir = args['datadir']
    infile = args['infile']
    mask = args['mask']
    labels = args['labels']
    defs = args['defs']
   

    # Import expression data -------------------------------------------------
    
    #File containing voxel expression matrix
    print("Importing gene-by-voxel expression matrix: {} ...".format(infile))
    
    #Import voxel expression matrix
    dfExprVoxel = pd.read_csv(os.path.join(datadir, infile))
    
    #Extract gene names from data frame
    genes = dfExprVoxel['Gene']

    #Extract voxels from data frame and convert to numpy array
    npExprVoxel = np.array(dfExprVoxel.loc[:,dfExprVoxel.columns != 'Gene'])
    
    #Transpose numpy array so that voxels are rows
    npExprVoxel = np.transpose(npExprVoxel)

    #Convert transposed numpy array to data frame.
    #Set genes to be column names.
    dfExprVoxelTranspose = pd.DataFrame(npExprVoxel, columns=genes)


    # Import imaging data ----------------------------------------------------
    
    #Import image mask, flatten and convert to numpy array
    print("Importing mask volume: {} ...".format(mask))

    maskVol = volumeFromFile(os.path.join('data', 'imaging', mask))
    maskArray = np.array(maskVol.data.flatten())
    maskVol.closeVolume()
    
    print("Importing DSURQE label volume: {} ...".format(labels))
    
    #Import DSURQE label volume, flatten, and convert to numpy array
    labelVol = volumeFromFile(os.path.join('data', 'imaging', labels))
    labelArray = np.array(labelVol.data.flatten())
    labelVol.closeVolume()

    #Mask the label array
    labelArrayMasked = labelArray[maskArray == 1]
    
    print("Importing DSURQE label definitions: {} ...".format(defs))
    
    #Import DSURQE label definitions
    dfAtlasDefs = pd.read_csv(os.path.join('data', 'imaging', defs))
    
    
    # Match labels to voxels -------------------------------------------------
        
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
    indNotMissing = dfExprVoxelTranspose['Region'] != 'nan'
    dfExprVoxelTranspose = dfExprVoxelTranspose.loc[indNotMissing]
    
    
    # Aggregate expression data ----------------------------------------------
    
    print("Aggregating expression data...")
    
    #Aggregate expression data by ROI
    dfExprRegionTranspose = (dfExprVoxelTranspose
                             .groupby('Region')
                             .aggregate(np.mean))

    #Transpose data frame so that rows are genes
    dfExprRegion = np.transpose(dfExprRegionTranspose)

    #Gene names are stored in the index
    #Label the index as Gene
    dfExprRegion.index.name = 'Gene'


    # Write ------------------------------------------------------------------

    print("Writing to file...")
    
    #Write regional expression matrix to csv
    dfExprRegion.to_csv(os.path.join(datadir, args['outfile']))
    
    return
    
if __name__ == '__main__':
    main()
