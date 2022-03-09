# build_voxel_matrix.py --------------------------------------------
#
#
#
# Antoine Beauchamp
# Created: June 19th, 2020
# Edited: March 7th, 2022
# --------------------------------------------------------------------

# Packages -----------------------------------------------------------

import argparse
import os
import re
import warnings
import numpy as np
import pandas as pd
import multiprocessing as mp

from pyminc.volumes.factory import *
from glob                   import glob
from tqdm                   import tqdm
from functools              import partial
from sklearn.impute         import KNNImputer
from sklearn.preprocessing  import FunctionTransformer
from sklearn.pipeline       import Pipeline


# Functions ----------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/expression/',
        help = 'Directory containing expression data'
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = 'Directory in which to save expression matrix'
    )
    
    parser.add_argument(
        '--dataset',
        type = str,
        default = 'coronal',
        choices = ['coronal', 'sagittal'],
        help = 'AMBA dataset to import'
    )
    
    parser.add_argument(
        '--mask',
        type = str,
        default = 'coronal',
        choices = ['coronal', 'sagittal'],
        help = 'Mask to apply to ISH images'
    )
    
    parser.add_argument(
        '--log2',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Option to transform expression data to log2'
    )
    
    parser.add_argument(
        '--groupexp',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Option to group multiple ISH experiments per gene'
    )
    
    parser.add_argument(
        '--threshold',
        type = float,
        default = 0.2,
        help = ''
    )
    
    parser.add_argument(
        '--impute',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Option to impute missing data using KNN imputation'
    )
    
    parser.add_argument(
        '--parallel',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Option to import image files in parallel'
    )
    
    parser.add_argument(
        '--nproc',
        type = int,
        default = mp.cpu_count(),
        help = 'Number of processors to use in parallel. Ignored if --parallel set to false.'
    )
    
    args = vars(parser.parse_args())
    
    return args
    
    
def importImage(img, mask):
    
    """ """
    
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


def buildExpressionMatrix(files, mask, log_transform = True, group_experiments = True, threshold = 0.2, parallel = True, nproc = None):
    
    """ """
    
    importImage_partial = partial(importImage, mask = mask)

    if parallel:

        if nproc is None:
            nproc = mp.cpu_count()

        pool = mp.Pool(nproc)

        arrays = []
        for array in tqdm(pool.imap(importImage_partial, files), total = len(files)):
            arrays.append(array)

        pool.close()
        pool.join()

    else:

        arrays = list(map(importImage_partial, tqdm(files)))
    
    dfExpression = pd.DataFrame(np.asarray(arrays), index = [os.path.basename(file) for file in files])

    #Transform to log2
    if log_transform:
        print("Applying log2 transform...")
        dfExpression = np.log2(dfExpression)
    
    dfExpression.index = dfExpression.index.str.replace('.mnc', '').str.replace('_.*', '')
        
    dfExpression.index.name = 'Gene'
    
    #Aggregate experiments per gene if flag is set
    if group_experiments:
        print("Aggregating multiple experiments per gene...")
        dfExpression = dfExpression.groupby(dfExpression.index).aggregate(np.mean)
    
    #Remove genes where a threshold of voxels aren't expressing
    fracVoxelsNA = dfExpression.isna().sum(axis=1)/len(dfExpression.columns)
    dfExpression = dfExpression[fracVoxelsNA < threshold]
    
    return dfExpression
    

def main():

    #Load command line arguments
    args = parse_args()
    datadir = args['datadir']
    outdir = args['outdir']
    dataset = args['dataset']
    mask = args['mask']
    
    print("Importing {} dataset using {} mask".format(dataset, mask))
    
    if (dataset == 'sagittal') and (mask == 'coronal'):
        warnings.warn("Running with sagittal dataset and coronal mask is not ideal. Proceed with caution.")
        
    
    #If dataset is sagittal, use only those genes that are also in the coronal set
    if dataset == "sagittal":
        
        #Paths to sagittal and coronal data set directories
        pathGeneDir_Sagittal = os.path.join(datadir, dataset, '')
        pathGeneDir_Coronal = os.path.join(datadir, 'coronal', '')
    
        #Build paths to all files in the directories
        pathGeneFiles_Sagittal = np.array(glob(pathGeneDir_Sagittal + "*.mnc"))
        pathGeneFiles_Coronal = np.array(glob(pathGeneDir_Coronal + "*.mnc"))
    
        #Extract gene names for coronal and sagittal data sets
        genes_Sagittal = np.array([re.sub(r"_[0-9]+.mnc", "", file) for file in 
                                    [os.path.basename(path) for path in pathGeneFiles_Sagittal]])
        genes_Coronal = np.array([re.sub(r"_[0-9]+.mnc", "", file) for file in 
                                    [os.path.basename(path) for path in pathGeneFiles_Coronal]])

        #Identify genes from sagittal data in coronal data
        isInCoronal = np.isin(genes_Sagittal, genes_Coronal)

        #Extract subset of sagittal gene files
        pathGeneFiles = pathGeneFiles_Sagittal[isInCoronal]
        
    else:
        pathGeneDir = os.path.join(datadir, dataset, '')
        pathGeneFiles = np.array(glob(pathGeneDir+"*.mnc"))

    print("Building voxel expression matrix...")

    #Mask files
    if mask == "sagittal":
        maskfile = "data/imaging/sagittal_200um_coverage_bin0.8.mnc"
    else: 
        maskfile = "data/imaging/coronal_200um_coverage_bin0.8.mnc"
    
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
        dfExpression = pd.DataFrame(imputing_pipeline.fit_transform(dfExpression.to_numpy()),
                                    index = genes)
        
        
        
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
