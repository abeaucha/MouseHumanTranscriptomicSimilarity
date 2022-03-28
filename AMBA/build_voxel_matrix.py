# ----------------------------------------------------------------------------
# build_voxel_matrix.py 
# Author: Antoine Beauchamp
# Created: June 19th, 2020

"""
Build a gene-by-voxel expression matrix

Description
-----------
This script imports in-situ hybridization (ISH) MINC images from the
Allen Mouse Brain Atlas, applies some pre-processing, and stores the
voxel-wise expression values for ISH experiments or genes in a matrix. 
This matrix is written out to a CSV file. 

Pre-processing options include:
1. log2 transformation
2. Voxel-wise average of the expression from multiple experiments
   that correspond to a single gene
3. Filtering out genes or experiments that have too many empty voxels
4. Imputing empty voxels using K-nearest neighbours

The script can import either the coronal or sagittal AMBA data sets,
using either a (bilateral) coronal or (unilateral) sagittal mask. When
importing the sagittal data set, the script will only import those
genes that are also present in the coronal data set.
"""

# Packages -------------------------------------------------------------------

import argparse
import os
import warnings
import numpy as np
import pandas as pd
import multiprocessing as mp

from pyminc.volumes.factory import *
from re                     import sub
from glob                   import glob
from tqdm                   import tqdm
from functools              import partial
from sklearn.impute         import KNNImputer
from sklearn.preprocessing  import FunctionTransformer
from sklearn.pipeline       import Pipeline

# Functions ------------------------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )

    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/expression/',
        help = ("Directory containing AMBA data. This directory should "
                "contain sub-directories 'coronal' and 'sagittal', which "
                "contain expression MINC files.")
    )
    
    parser.add_argument(
        '--imgdir',
        type = str,
        default = 'data/imaging/',
        help = ("Directory containing imaging data.")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = "Directory in which to save the expression matrix."
    )
    
    parser.add_argument(
        '--dataset',
        type = str,
        default = 'coronal',
        choices = ['coronal', 'sagittal'],
        help = "AMBA dataset to import."
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        default = 'coronal',
        choices = ['coronal', 'sagittal'],
        help = "Mask to apply to expression images."

    )
    
    parser.add_argument(
        '--log2',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = "Option to transform expression data to log2."
    )
    
    parser.add_argument(
        '--groupexp',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = "Option to group multiple ISH experiments per gene."
    )
    
    parser.add_argument(
        '--threshold',
        type = float,
        default = 0.2,
        help = ("Maximal fraction of empty voxels allowed in an ISH "
                "image to be included in the final matrix.")
    )
    
    parser.add_argument(
        '--impute',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Option to impute empty voxels using "
                "K-nearest neighbours imputation.")
    )
    
    parser.add_argument(
        '--parallel',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = "Option to import image files in parallel."
    )
    
    parser.add_argument(
        '--nproc',
        type = int,
        default = mp.cpu_count(),
        help = ("Number of processors to use in parallel. "
                "Ignored if --parallel set to false.")
    )
    
    args = vars(parser.parse_args())
    
    return args
    
    
def importImage(img, mask):
    
    """
    Import a MINC file as NumPy array

    Description
    -----------
    This function imports a MINC file as a flattened NumPy array. 
    The image is masked using the mask provided, and values of -1
    and 0 are replaced with NumPy NaNs.

    Arguments
    ---------
    img: str
        Path to the MINC file to import.
    mask: str
        Path to the the MINC file containing the mask. Must be in
        the same space as `img`.

    Returns
    -------
    imageArrayMasked: numpy.ndarray
        A 1-dimensional NumPy array containing the masked image voxel
        values.
    """
    
    #Import mask mask and convert to numpy array
    maskVol = volumeFromFile(mask)
    maskArray = np.array(maskVol.data.flatten())
    maskVol.closeVolume()
    
    #Read ISH data to numpy array
    imageVol = volumeFromFile(img)
    imageArray = np.array(imageVol.data.flatten())
    imageVol.closeVolume()
    
    #Apply mask and convert -1 to NaN
    imageArrayMasked = imageArray[maskArray == 1]
    imageArrayMasked[imageArrayMasked == -1] = np.nan
    imageArrayMasked[imageArrayMasked == 0] = np.nan
    
    return imageArrayMasked


def buildExpressionMatrix(files, mask, log_transform = True,
                          group_experiments = True, threshold = 0.2, 
                          parallel = True, nproc = None):
    
    """ 
    Build gene-by-voxel expression matrix

    Arguments
    ---------
    files: list of str
        List containing paths to expression MINC files.
    mask: str
        Path to mask MINC file.
    log_transform: bool, optional
        Option to apply a log2 transform to the expression values.
        (default True)
    group_experiments: bool, optional,
        Option to compute the voxel-wise average of expression values 
        for experiments that correspond to the same gene. (default True)
    threshold: float, optional
        Threshold value indicating the fraction of empty voxels in an 
        image above which the image is discarded (default 0.2)
    parallel: bool, optional
        Option to import MINC files in parallel. (default True)
    nproc: int, optional
        Number of CPUs to use in parallel. If `None`, all CPUs are
        used. (default None)

    Returns
    -------
    dfExpression: pandas.core.frame.DataFrame
       A DataFrame containing the expression of experiments/genes in
       the Allen Mouse Brain Atlas. If `group_experiments` is False,
       every row corresponds to an experiment. If `group_experiments`
       is True, every row corresponds to a gene.
    """
    
    importImage_partial = partial(importImage, mask = mask)

    if parallel:

        if nproc is None:
            nproc = mp.cpu_count()

        pool = mp.Pool(nproc)

        arrays = []
        for array in tqdm(pool.imap(importImage_partial, files),
                          total = len(files)):
            arrays.append(array)

        pool.close()
        pool.join()

    else:

        arrays = list(map(importImage_partial, tqdm(files)))
    
    dfExpression = pd.DataFrame(np.asarray(arrays), 
                   index = [os.path.basename(file) for file in files])

    #Transform to log2
    if log_transform:
        print("Applying log2 transform...")
        dfExpression = np.log2(dfExpression)
    
    dfExpression.index = (dfExpression.index
                          .str.replace('.mnc', '')
                          .str.replace('_.*', ''))
        
    dfExpression.index.name = 'Gene'
    
    #Aggregate experiments per gene if flag is set
    if group_experiments:
        print("Aggregating multiple experiments per gene...")
        dfExpression = (dfExpression
                        .groupby(dfExpression.index)
                        .aggregate(np.mean))
    
    #Remove genes where a threshold of voxels aren't expressing
    fracVoxelsNA = dfExpression.isna().sum(axis=1)/len(dfExpression.columns)
    dfExpression = dfExpression[fracVoxelsNA < threshold]
    
    return dfExpression

# Main -----------------------------------------------------------------------

def main():

    #Load command line arguments
    args = parse_args()
    datadir = args['datadir']
    imgdir = args['imgdir']
    outdir = args['outdir']
    dataset = args['dataset']
    mask = args['mask']
    
    print("Importing {} dataset using {} mask".format(dataset, mask))
    
    if (dataset == 'sagittal') and (mask == 'coronal'):
        warnings.warn(("Running with sagittal dataset and coronal mask is "
                       "not ideal. Proceed with caution."))
    
    #If dataset is sagittal, use only those genes that are also in the
    #coronal set
    if dataset == "sagittal":
        
        #Paths to sagittal and coronal data set directories
        pathGeneDir_Sagittal = os.path.join(datadir, dataset, '')
        pathGeneDir_Coronal = os.path.join(datadir, 'coronal', '')

        #Build paths to all files in the directories
        pathGeneFiles_Sagittal = glob(pathGeneDir_Sagittal + '*.mnc')
        pathGeneFiles_Coronal = glob(pathGeneDir_Coronal + '*.mnc')

        #Extract gene names for coronal and sagittal data sets
        genes_Sagittal = [sub(r'_[0-9]+.mnc', '', file) for file in 
                [os.path.basename(path) for path in pathGeneFiles_Sagittal]]
        genes_Coronal = [sub(r'_[0-9]+.mnc', "", file) for file in 
                [os.path.basename(path) for path in pathGeneFiles_Coronal]]

        #Identify genes from sagittal data in coronal data
        isInCoronal = np.isin(np.array(genes_Sagittal),
                              np.array(genes_Coronal))

        #Extract subset of sagittal gene files
        pathGeneFiles_Sagittal = np.array(pathGeneFiles_Sagittal)
        pathGeneFiles = list(pathGeneFiles_Sagittal[isInCoronal])
        
    else:
        pathGeneDir = os.path.join(datadir, dataset, '')
        pathGeneFiles = glob(pathGeneDir+'*.mnc')

    print("Building voxel expression matrix...")

    #Mask files
    if mask == 'sagittal':
        maskfile = os.path.join(imgdir, 'sagittal_200um_coverage_bin0.8.mnc')
    else: 
        maskfile = os.path.join(imgdir, 'coronal_200um_coverage_bin0.8.mnc')
    
    #Build expression data frame
    log_transform = True if args['log2'] == 'true' else False
    groupexp = True if args['groupexp'] == 'true' else False
    parallel = True if args['parallel'] == 'true' else False
    dfExpression = buildExpressionMatrix(files = pathGeneFiles, 
                                         mask = maskfile,
                                         log_transform = log_transform,
                                         group_experiments = groupexp, 
                                         threshold = args['threshold'], 
                                         parallel = parallel, 
                                         nproc = args['nproc'])

    #Impute missing values
    impute = True if args['impute'] == 'true' else False
    if impute:
        
        print("Imputing missing values using K-nearest neighbours...")
        
        #Initialize imputer and transposer
        imputer = KNNImputer(missing_values = np.nan)
        transposer = FunctionTransformer(np.transpose)
        
        #Build pipeline
        imputing_pipeline = Pipeline([('transpose1', transposer),
                                     ('impute', imputer),
                                     ('transpose2', transposer)])
        
        #Store gene names
        genes = dfExpression.index
        
        #Impute missing values and assign as data frame
        dfExpression = pd.DataFrame(imputing_pipeline.fit_transform(
                                                  dfExpression.to_numpy()
                                    ), index = genes)

    #Write to file
    print("Writing to file...")
    
    outfile = 'MouseExpressionMatrix_voxel_{}_mask{}'.format(dataset, mask)
    
    if log_transform:
        outfile = outfile+'_log2'
        
    if groupexp:
        outfile = outfile+'_grouped'
        
    if impute:
        outfile = outfile+'_imputed'
    
    outfile = outfile+'.csv'
    
    dfExpression.to_csv(os.path.join(outdir, outfile))

    return

if __name__ == '__main__':
    main()
