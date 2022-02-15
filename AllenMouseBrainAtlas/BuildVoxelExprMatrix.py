# BuildVoxelExprMatrix.py --------------------------------------------
#
#
#
# Antoine Beauchamp
# Created: June 19th, 2020
# Edited: February 15th, 2022
# --------------------------------------------------------------------

# Packages -----------------------------------------------------------

import numpy as np
import pandas as pd
import re
import argparse
import sys
import warnings

from pyminc.volumes.factory import *
from glob                   import glob
from os.path                import basename

from sklearn.impute         import KNNImputer
from sklearn.preprocessing  import FunctionTransformer
from sklearn.pipeline       import Pipeline


# Functions ----------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser()
    
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
        '--groupexp',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Option to group multiple ISH experiments per gene'
    )
    
    parser.add_argument(
        '--impute',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Option to impute missing data using KNN imputation'
    )
    
    args = vars(parser.parse_args())
    
    return args
    
    

def main():

    #Load command line arguments
    args = parse_args()
    
    print('Running with options: {}'.format(args))
    
    #Flags
    dataset = args['dataset']
    maskFlag = args['mask']
    groupExperiments = True if args['groupexp'] == 'true' else False
    impute = True if args['impute'] == 'true' else False

    
    print("Using {} dataset with {} mask".format(dataset, maskFlag))
    
    if (dataset == 'sagittal') and (maskFlag == 'coronal'):
        warnings.warn("Running with sagittal dataset and coronal mask is not ideal. Proceed with caution.")
        
    
    #Paths
    pathAMBA = "/projects/abeauchamp/Projects/MouseHumanMapping/Paper_TranscriptomicSimilarity/AllenMouseBrainAtlas/
    
    #If dataset is sagittal, use only those genes that are also in the coronal set
    if dataset == "sagittal":
        
        #Paths to sagittal and coronal data set directories
        pathGeneDir_Sagittal = "/projects/yyee/tools/Allen_brain/data/expression/P56/sagittal/"
        pathGeneDir_Coronal = "/projects/yyee/tools/Allen_brain/data/expression/P56/coronal/"
    
        #Build paths to all files in the directories
        pathGeneFiles_Sagittal = np.array(glob(pathGeneDir_Sagittal + "*.mnc"))
        pathGeneFiles_Coronal = np.array(glob(pathGeneDir_Coronal + "*.mnc"))
    
        #Extract gene names for coronal and sagittal data sets
        genes_Sagittal = np.array([re.sub(r"_[0-9]+.mnc", "", file) for file in 
                                    [basename(path) for path in pathGeneFiles_Sagittal]])
        genes_Coronal = np.array([re.sub(r"_[0-9]+.mnc", "", file) for file in 
                                    [basename(path) for path in pathGeneFiles_Coronal]])

        #Identify genes from sagittal data in coronal data
        isInCoronal = np.isin(genes_Sagittal, genes_Coronal)

        #Extract subset of sagittal gene files
        pathGeneFiles = pathGeneFiles_Sagittal[isInCoronal]
        
    else:
        pathGeneDir = "/projects/yyee/tools/Allen_brain/data/expression/P56/coronal/"
        pathGeneFiles = np.array(glob(pathGeneDir + "*.mnc"))


    #Load image mask and convert to numpy array
    if maskFlag == "sagittal":
        maskVol = volumeFromFile(pathAMBA+'Data/MRI/sagittal_200um_coverage_bin0.8.mnc')        
    else: 
        maskVol = volumeFromFile(pathAMBA+'Data/MRI/coronal_200um_coverage_bin0.8.mnc')
        
    maskArray = np.array(maskVol.data.flatten())
    maskVol.closeVolume()


    print("Building expression matrix...")
    
    #Build expression data frame
    dfExpression = buildVoxelExprMatrix(pathGeneFiles, mask = maskArray, groupExperiments = groupExperiments)

    
    #Switch to impute missing values
    if impute == True:
        
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
        
        imputed = '_imputed'
        
    else:
        imputed = ''
        
        
    
    print("Writing to file...")
    
    #Write to file
    if groupExperiments == False:
        wDups = '_wDups'
    else:
        wDups = ''
    
    outfile = 'MouseExpressionMatrix_Voxel_'+dataset+'_mask'+maskFlag+wDups+imputed+'.csv'
    
    
    dfExpression.to_csv(pathAMBA+'Data/'+outfile)
    
    

def buildVoxelExprMatrix(paths, mask, threshold = 0.2, groupExperiments = True):
    
    """ """
    
    #Load expression data
    dfExpression = loadExpressionData(paths, mask)

    #Transform to log2 and normalize by experiment
    dfExpression = np.log2(dfExpression)
    
    #Replace index with gene acronyms
    dfExpression.index = dfExpression.index.str.replace('.mnc', '').str.replace('_.*', '')
    dfExpression.index.name = 'Gene'
    
    #Aggregate experiments per gene if flag is set
    if groupExperiments == True:
        print("Aggregating multiple experiments per gene...")
        dfExpression = dfExpression.groupby(dfExpression.index).aggregate(np.mean)
    
    #Remove genes where a threshold of voxels aren't expressing
    fracVoxelsNA = dfExpression.isna().sum(axis=1)/len(dfExpression.columns)
    dfExpression = dfExpression[fracVoxelsNA < threshold]
    
    return dfExpression
    

    
def loadExpressionData(paths, mask):

    """ """

    #Number of files and number of voxels per file
    nFiles = len(paths)
    nVoxels = int(np.sum(mask))

    #Initialize matrix
    matExpr = np.empty((nFiles, nVoxels), dtype = 'float')

    #Iterate over files
    for i, path in enumerate(paths):

        if (i%100 == 0):
            print('On file {} of {}'.format(i, len(paths)))

        #Read ISH data to numpy array
        exprVol = volumeFromFile(path)
        exprArray = np.array(exprVol.data.flatten())
        exprVol.closeVolume()

        #Apply mask and convert -1 to NaN
        exprArrayMasked = exprArray[mask == 1]
        exprArrayMasked[exprArrayMasked == -1] = np.nan
        exprArrayMasked[exprArrayMasked == 0] = np.nan

        #Write expression data to matrix
        matExpr[i,] = exprArrayMasked
        
    #Convert matrix to df, using file names as index
    dfExpression = pd.DataFrame(matExpr, index = [basename(path) for path in paths])

    return dfExpression


if __name__ == '__main__':
    main()
