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
import numpy                as np
import pandas               as pd
from pyminc.volumes.factory import volumeFromFile

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
        '--imgdir',
        type = str,
        default = 'data/imaging/',
        help = "Directory containing imaging files"
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
                "voxel-wise expression matrix provided to --infile. "
                "Must reside in --imgdir.")
    )
    
    parser.add_argument(
        '--labels',
        type = str,
        help = ("Name of MINC file containing the atlas labels to use. "
                "Must be in CCFv3 space. Must reside in --imgdir.")
    )
    
    parser.add_argument(
        '--defs',
        type = str,
        help = ("Name of CSV file containing the names of the "
                "neuroanatomical regions corresponding to the "
                "atlas labels in --labels. Must reside in --imgdir.")
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Verbosity.'
    )
    
    args = vars(parser.parse_args())
    
    return args

# Main -----------------------------------------------------------------------

def main():
    
    #Load command line arguments
    args = parse_args()
    
    #Flags
    datadir = args['datadir']
    imgdir = args['imgdir']
    infile = args['infile']
    outfile = args['outfile']
    mask = args['mask']
    labels = args['labels']
    defs = args['defs']
    
    if infile is None:
        raise Exception("No input file passed to argument --infile")

    if outfile is None:
        raise Exception("No output file name passed to argument --outfile")
        
    if mask is None:
        raise Exception("No mask file passed to argument --mask")
        
    if labels is None:
        raise Exception("No label file passed to argument --labels")
        
    if defs is None:
        raise Exception("No atlas definitions file passed to argument --defs")
        
    datadir = os.path.join(datadir, '')
    
 
    # Import expression data -------------------------------------------------
    
    #File containing voxel expression matrix
    
    try: 
        if verbose:
            print("Importing gene-by-voxel expression matrix: {} ...".format(infile))
    
        #Import voxel expression matrix
        dfExprVoxel = pd.read_csv(os.path.join(datadir, infile))
    except FileNotFoundError:
        raise FileNotFoundError("Input file {} not found in data directory {}"
                               .format(infile, datadir))
        
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
    
    imgdir = os.path.join(imgdir, '')
    
    if os.path.exists(os.path.join(imgdir, mask)) == False:
        raise FileNotFoundError("Mask file {} not found in imaging directory {}"
                               .format(mask, imgdir))
    
    #Import image mask, flatten and convert to numpy array
    if verbose:
        print("Importing mask volume: {} ...".format(mask))

    maskVol = volumeFromFile(os.path.join(imgdir, mask))
    maskArray = np.array(maskVol.data.flatten())
    maskVol.closeVolume()
    
    if os.path.exists(os.path.join(imgdir, labels)) == False:
        raise FileNotFoundError("Label file {} not found in imaging directory {}"
                               .format(labels, imgdir))
    
    if verbose:
        print("Importing atlas label volume: {} ...".format(labels))
    
    #Import DSURQE label volume, flatten, and convert to numpy array
    labelVol = volumeFromFile(os.path.join(imgdir, labels))
    labelArray = np.array(labelVol.data.flatten())
    labelVol.closeVolume()

    #Mask the label array
    labelArrayMasked = labelArray[maskArray == 1]
    
    if os.path.exists(os.path.join(imgdir, defs)) == False:
        raise FileNotFoundError("Atlas definitions file {} not found in imaging "
                                "directory {}".format(labels, imgdir)) 
    
    if verbose:
        print("Importing atlas label definitions: {} ...".format(defs))
    
    #Import atlas label definitions
    dfAtlasDefs = pd.read_csv(os.path.join(imgdir, defs))
    
    
    # Match labels to voxels -------------------------------------------------
        
    if verbose:
        print("Matching atlas labels to voxels...")

    #Initialize empty string array to match atlas regions to voxels
    structArray = np.empty(len(labelArrayMasked), dtype = 'U100')

    #Iterate over atlas regions and match to voxels
    for i, row in dfAtlasDefs.iterrows():
        lab = row['Label']
        indlab = labelArrayMasked == row['Label']
        
        if lab == 0:
            structArray[indlab] = np.nan
        else:
            structArray[indlab] = row['Structure']
    
    #Assign atlas regions to new column in the data frame
    dfExprVoxelTranspose['Region'] = structArray
    
    #Remove voxels that are empty
    indNotMissing = dfExprVoxelTranspose['Region'] != 'nan'
    dfExprVoxelTranspose = dfExprVoxelTranspose.loc[indNotMissing]
    
    
    # Aggregate expression data ----------------------------------------------
    
    if verbose:
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

    if verbose:
        print("Writing to file...")

    #Write regional expression matrix to csv
    dfExprRegion.to_csv(os.path.join(datadir, outfile))
    
    return
    
if __name__ == '__main__':
    main()
