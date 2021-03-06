# ----------------------------------------------------------------------------
# Model_MultilayerPerceptron_Validation_CoronalSagittalSampling.py
# Author: Antoine Beauchamp
# Created: February 3rd, 2021

"""
Implement cross-validation 

Description
-----------
"""

# Packages -------------------------------------------------------------------

import pandas as pd
import numpy as np
import pickle
import argparse
import sys

from itertools                import product

from sklearn.impute           import SimpleImputer
from sklearn.preprocessing    import StandardScaler, FunctionTransformer
from sklearn.pipeline         import Pipeline
from sklearn.metrics          import accuracy_score, confusion_matrix

from skorch                   import NeuralNetClassifier
from skorch.toy               import make_classifier
from skorch.helper            import DataFrameTransformer
from skorch.callbacks         import LRScheduler
from torch.optim              import AdamW
from torch.optim.lr_scheduler import OneCycleLR
from torch.cuda               import is_available


# Functions ------------------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--labels",
        type = str,
        default = 'region5',
        choices = ['region5', 
                   'region11', 
                   'region28', 
                   'region46', 
                   'region67', 
                   'region130'],
        help = "Class of labels on which to train."
    )
    
    parser.add_argument(
        "--nunits",
        nargs = "*",
        type = int,
        default = [500],
        help = "List containing the number of hidden units to tune over."
    )

    parser.add_argument(
        "--nlayers",
        nargs = "*",
        type = int,
        default = [4],
        help = "List containing the number of hidden layers to tune over."
    )

    parser.add_argument(
        "--dropout",
        nargs = "*",
        type = float,
        default = [0.],
        help = "List containing dropout rates to tune over."
    )

    parser.add_argument(
        "--L2",
        nargs = "*",
        type = float,
        default = [1e-6],
        help = "List containing weight decay values to tune over."
    )

    parser.add_argument(
        "--nsamples",
        type = int,
        default = 1,
        help = "Number of times to train and evaluate each hyperparameter combination."
    )

    parser.add_argument(
        "--nepochs",
        type = int,
        default = 200,
        help = "Number of epochs to train over."
    )

    parser.add_argument(
        "--learningrate",
        type = float,
        default = 1e-5,
        help = "Learning rate during training."
    )

    parser.add_argument(
        "--confusionmatrix",
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = "Flag to indicate whether to compute confusion matrices."
    )
    
    parser.add_argument(
        "--imputed",
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = "Flag to indicate whether to use imputed data or data with missing values. If false, will impute missing values using median imputation [default true]"
    )
    
    args = vars(parser.parse_args())

    return args
    
    
def buildTrainValidationSets(coronal, sagittal, seed = None):
    
    """ """
    
    #Get genes in the coronal data set (includes duplicates)
    genes_coronal = coronal.columns.str.replace('\.\.\.[0-9]+', '')
    
    #Get unique genes, present in sagittal and coronal data
    genes_unique = np.unique(genes_coronal)
    
    #Initialize data frames for training and validation sets
    train = pd.DataFrame(np.empty((coronal.shape[0], len(genes_unique)), dtype = 'float'),
                         columns = genes_unique)
    
    validation = pd.DataFrame(np.empty_like(train),
                              columns = genes_unique)
    
    #Initialize random number generator
    rng = np.random.default_rng(seed = seed)
    
    #Iterate over unique genes
    for i, gene in enumerate(genes_unique):
        
        #If the gene has a replicated coronal experiment, use that
        if np.sum(genes_coronal == gene) > 1:
            
            #Extract replicated experiments
            df_choices = coronal.loc[:, genes_coronal == gene]
            
            #Randomly choose one of the experiments for the training set
            choices = np.arange(0, df_choices.shape[1])
            choice_train = rng.choice(choices, 1)
            
            #Randomly choose one of the remaining experiments for the validation set
            choice_val = rng.choice(choices[choices != choice_train], 1)
            
            #Assign data to training and validation sets
            train.iloc[:,i] = df_choices.iloc[:, choice_train]
            validation.iloc[:, i] = df_choices.iloc[:, choice_val]
            
        #If the gene is unique in the coronal set, choose between coronal and sagittal
        else:
            
            #Random binary choice
            choice_train = rng.choice([0,1],1)
            
            #Assign coronal data to training and sagittal to validation, or vice versa
            if choice_train[0] == 0:
                train.iloc[:,i] = coronal.loc[:, coronal.columns == gene]
                validation.iloc[:,i] = sagittal.loc[:, sagittal.columns == gene]
            else:
                train.iloc[:,i] = sagittal.loc[:, sagittal.columns == gene]
                validation.iloc[:,i] = coronal.loc[:, coronal.columns == gene]
                
    return train, validation


def main():

    #Load command line arguments
    args = parse_args()
   
    print('Running with options: {}'.format(args))
    
    
    # Importing data -------------------------------------------------------------------

    #Paths
    pathHome = "/hpf/largeprojects/MICe/abeauchamp/Projects/MouseHumanMapping/Paper_Descriptive/"

    #Set up filepaths
    if args['imputed'] == 'true':
        fileSagittal = "MouseExpressionMatrix_Voxel_sagittal_masksagittal_imputed_labelled.csv"
        fileCoronal = "MouseExpressionMatrix_Voxel_coronal_masksagittal_wDups_imputed_labelled.csv"
    else: 
        fileSagittal = "MouseExpressionMatrix_Voxel_sagittal_masksagittal_labelled.csv"
        fileCoronal = "MouseExpressionMatrix_Voxel_coronal_masksagittal_wDups_labelled.csv"

    pathSagittal = pathHome+"Data/"+fileSagittal
    pathCoronal = pathHome+"Data/"+fileCoronal

    print("Importing coronal and sagittal data sets...")

    #Import data
    dfSagittal = pd.read_csv(pathSagittal)
    dfCoronal = pd.read_csv(pathCoronal)



    # Processing -------------------------------------------------------------------------

    print("Cleaning and preparing data...")

    #Identify and remove label columns
    indLabelsCoronal = dfCoronal.columns.str.match('Region')
    dfInputCoronal = dfCoronal.loc[:, ~indLabelsCoronal]

    indLabelsSagittal = dfSagittal.columns.str.match('Region')
    dfInputSagittal = dfSagittal.loc[:, ~indLabelsSagittal]

    #Identify and remove columns in the sagittal dataset that contain infinities
    if args['imputed'] == 'false':
        whichInf = []
        for i in range(dfInputSagittal.shape[1]):
            ninf = np.sum(np.isinf(dfInputSagittal.iloc[:,i]))
            if(ninf > 0):
                whichInf.append(i)

        colsInf = dfInputSagittal.columns[whichInf]
        dfInputSagittal = dfInputSagittal.loc[:, ~dfInputSagittal.columns.isin(colsInf)]


    #Get genes in coronal and sagittal sets
    genesCoronal = dfInputCoronal.columns.str.replace('\.\.\.[0-9]+', '')
    genesSagittal = dfInputSagittal.columns.str.replace('\.\.\.[0-9]+', '')

    #Extract genes in the sagittal set that are in the coronal set
    indSagittalInCoronal = genesSagittal.isin(genesCoronal)
    dfInputSagittal = dfInputSagittal.loc[:, indSagittalInCoronal]

    #Get new set of genes for sagittal set
    genesSagittal = dfInputSagittal.columns.str.replace('\.\.\.[0-9]+', '')

    #Extract genes in the coronal set that are in the sagittal set
    indCoronalInSagittal = genesCoronal.isin(genesSagittal)
    dfInputCoronal = dfInputCoronal.loc[:, indCoronalInSagittal]

    #Get new set of genes for coronal set
    genesCoronal = dfInputCoronal.columns.str.replace('\.\.\.[0-9]+', '')

        
    #Extract labels
    labelcol = args['labels'].title()
    dfLabels = dfCoronal[[labelcol]].copy()

    #Convert labels to category
    dfLabels.loc[:,labelcol] = dfLabels.loc[:,labelcol].astype('category')

    
    
    # ---------

    print("Beginning training and validation...")

    #Define a dictionary containing the grid values
    dictGrid = {'sample':np.arange(0,args['nsamples']),
               'hidden_units':args['nunits'],
               'hidden_layers':args['nlayers'],
               'dropout':args['dropout'],
               'weight_decay':args['L2']}
    

    #Expand the dictionary grid into a data frame containing all combinations
    dfPerformance = pd.DataFrame([row for row in product(*dictGrid.values())], columns = dictGrid.keys())

    #Initialize training and validation accuracy values
    dfPerformance['train_accuracy'] = 0
    dfPerformance['val_accuracy'] = 0

    #Set max number of epochs to train, and learning rate
    max_epochs = args['nepochs']
    learning_rate = args['learningrate']

    #Iterate over unique samples
    for sample in np.unique(dfPerformance['sample']):

        print('On sample {}'.format(sample))

        print('Generating training and validation sets...')

        #For the given sample, build the training and validation sets
        dfTraining, dfValidation = buildTrainValidationSets(dfInputCoronal, dfInputSagittal, seed = sample)    


        print('Preprocessing data...')

        #Initialize classes for imputing, scaling, centering and transposing
        scale = StandardScaler()
        center = StandardScaler(with_std = False)
        transpose = FunctionTransformer(np.transpose)
        
        if args['imputed'] == 'false':
            impute = SimpleImputer(missing_values = np.nan, strategy = 'median')
            processing_pipeline = Pipeline([('impute', impute),
                                            ('transpose1', transpose),
                                            ('scale', scale),
                                            ('transpose2', transpose),
                                            ('center', center)])
        else:
            processing_pipeline = Pipeline([('transpose1', transpose),
                                            ('scale', scale),
                                            ('transpose2', transpose),
                                            ('center', center)])

        #Fit the pipeline to the training data and transform
        X_train = processing_pipeline.fit_transform(dfTraining.to_numpy())
        X_train = X_train.astype(np.float32)

        #Transform labels into data preferred by the network
        y_train = DataFrameTransformer().fit_transform(dfLabels)[labelcol]

        #For the given training/validation sample, extract hyperparameters to iterate over
        dfParams = dfPerformance[dfPerformance['sample'] == sample]

        #Iterate over hyperparameter combinations
        for index, row in dfParams.iterrows():

            print('Index {}'.format(index))

            #Extract hyperparameters
            hidden_units = int(row['hidden_units'])
            hidden_layers = int(row['hidden_layers'])
            dropout = row['dropout']
            weight_decay = row['weight_decay']

            print('Labels: {}; Hidden units: {}; Hidden layers: {}; Dropout: {}; L2: {}'.format(args['labels'], hidden_units, hidden_layers, dropout, weight_decay))


            #Generate classifier module with specified architecture
            MLPModule = make_classifier(input_units = X_train.shape[1],
                                    output_units = len(np.unique(y_train)),
                                    hidden_units = hidden_units,
                                    num_hidden = hidden_layers,
                                    dropout = dropout)

            np.random.seed(123)

            net = NeuralNetClassifier(
                            MLPModule,
                            train_split = None,
                            optimizer = AdamW,
                            optimizer__weight_decay = weight_decay, #Weight decay regularization
                            max_epochs = max_epochs,
                            callbacks = [('lr_scheduler',
                                          LRScheduler(policy=OneCycleLR,
                                                      total_steps=max_epochs,  # Total number of iterations
                                                      cycle_momentum=False,    # AdamW doesn't support
                                                      max_lr=learning_rate))]  # Max learning rate
                        )



            if is_available() == True:
                print("Training network using GPU...")
            else:
                print("Training network using CPU...")

            #Fit the network to the training data
            net.fit(X_train, y_train)

            #Predict training labels
            y_train_pred = net.predict(X_train)

            #Compute and save training accuracy
            dfPerformance.loc[index, 'train_accuracy'] = accuracy_score(y_train, y_train_pred)

            print("Computing validation score...")

            #Process the validation data using the preprocessing pipeline
            X_val = processing_pipeline.transform(dfValidation.to_numpy())
            X_val = X_val.astype(np.float32)

            #Predict validation labels
            y_val_pred = net.predict(X_val)

            #Compute and save validation accuracy
            dfPerformance.loc[index, 'val_accuracy'] = accuracy_score(y_train, y_val_pred)
            
            if args['confusionmatrix'] == 'true':
                
                print("Computing confusion matrices...")
                
                dfConfusionMat_Train = pd.DataFrame(confusion_matrix(y_train, y_train_pred))
                dfConfusionMat_Val = pd.DataFrame(confusion_matrix(y_train, y_val_pred))
                
                dfLabelsDummy = dfLabels.copy()
                dfLabelsDummy['DummyVariable'] = y_train
                dfLabelsUnique = dfLabelsDummy.sort_values('DummyVariable').drop_duplicates()
                
                dfConfusionMat_Train.columns = dfLabelsUnique[labelcol].astype('str')
                dfConfusionMat_Train['TrueLabels'] = dfLabelsUnique[labelcol].astype('str').reset_index(drop = True)
                
                dfConfusionMat_Val.columns = dfLabelsUnique[labelcol].astype('str')
                dfConfusionMat_Val['TrueLabels'] = dfLabelsUnique[labelcol].astype('str').reset_index(drop = True)

                if args['imputed'] == 'false':
                    fileTrain = "MLP_Validation_CoronalSagittalSampling_ConfusionMatrix_Training_"+\
                    args['labels'].title()+\
                    "_ImputeMedians"+\
                    "_Sample"+str(sample)+\
                    "_Layers"+str(hidden_layers)+\
                    "_Units"+str(hidden_units)+\
                    "_Dropout"+str(dropout)+\
                    "_L2"+str(weight_decay)+".csv"

                    fileVal = "MLP_Validation_CoronalSagittalSampling_ConfusionMatrix_Validation_"+\
                    args['labels'].title()+\
                    "_ImputeMedians"+\
                    "_Sample"+str(sample)+\
                    "_Layers"+str(hidden_layers)+\
                    "_Units"+str(hidden_units)+\
                    "_Dropout"+str(dropout)+\
                    "_L2"+str(weight_decay)+".csv"
                
                else:
                    
                    fileTrain = "MLP_Validation_CoronalSagittalSampling_ConfusionMatrix_Training_"+\
                    args['labels'].title()+\
                    "_ImputeKNN"+\
                    "_Sample"+str(sample)+\
                    "_Layers"+str(hidden_layers)+\
                    "_Units"+str(hidden_units)+\
                    "_Dropout"+str(dropout)+\
                    "_L2"+str(weight_decay)+".csv"

                    fileVal = "MLP_Validation_CoronalSagittalSampling_ConfusionMatrix_Validation_"+\
                    args['labels'].title()+\
                    "_ImputeKNN"+\
                    "_Sample"+str(sample)+\
                    "_Layers"+str(hidden_layers)+\
                    "_Units"+str(hidden_units)+\
                    "_Dropout"+str(dropout)+\
                    "_L2"+str(weight_decay)+".csv"

                
                
                dfConfusionMat_Train.to_csv(pathHome+"Data/ValidationExperiments/"+fileTrain, 
                                            index = False)
                dfConfusionMat_Val.to_csv(pathHome+"Data/ValidationExperiments/"+fileVal, 
                                          index = False)
                
            
    
    if args['imputed'] == 'false':
        fileout = "MLP_Validation_CoronalSagittalSampling_"+\
                  args['labels'].title()+\
                  "_ImputeMedians_temp.csv"
    else:
        fileout = "MLP_Validation_CoronalSagittalSampling_"+\
                  args['labels'].title()+\
                  "_ImputeKNN_temp.csv"
        
    dfPerformance.to_csv(pathHome+"Data/ValidationExperiments/"+fileout, index=False)
    
    
if __name__ == "__main__":
    main()
